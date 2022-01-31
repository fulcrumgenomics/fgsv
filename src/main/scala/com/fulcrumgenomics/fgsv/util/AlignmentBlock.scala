package com.fulcrumgenomics.fgsv.util

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.util.LazyLogging
import enumeratum.EnumEntry

import scala.annotation.tailrec
import scala.collection.immutable.IndexedSeq
import scala.collection.{BitSet, mutable}

/** Represents an alignment block of read to the reference. The alignment may not consume all the read bases.
 *
 * @param origin         from which read (or both) this block comes from
 * @param start          the 1-based start in the read **in sequencing order**
 * @param end            the 1-based inclusive end in the read **in sequencing order**
 * @param positiveStrand true if the read was mapped to the positive strand, false otherwise
 * @param cigar          the cigar **in sequencing order**
 * @param range          the genomic range to which this alignment maps
 */
case class AlignmentBlock(origin: BlockOrigin, start: Int, end: Int, positiveStrand: Boolean, cigar: Cigar, range: GenomicRange) {
  require(0 < start)
  require(start <= end)

  /** Returns a bit set where each bit is a query base (in sequencing order), and bit is set if this alignment maps a base */
  def toQueryBitSet: BitSet = {
    scala.collection.BitSet(Range.inclusive(start=this.start, end=this.end): _*)
  }

  /** Returns if this alignment overlaps the other alignment on the genome, and maps to the same strand. */
  def overlapsWithStrand(other: AlignmentBlock): Boolean = {
    this.range.overlaps(other.range) && this.positiveStrand == other.positiveStrand
  }

  /** Merges the other block with itself, returning a new alignment block.  The two blocks must overlap on the genome
   * and map to teh same strand.  If merging blocks from different reads (i.e. R1 and R2), the block origin will be set
   * to both.  The resulting cigar will be empty, and the start and end offset in the read will both be 1.
   * */
  def merge(other: AlignmentBlock): AlignmentBlock = {
    require(this.overlapsWithStrand(other))
    val range = this.range.union(other.range)
    val origin = if (this.origin == other.origin) this.origin else BlockOrigin.Both
    AlignmentBlock(origin=origin, start=1, end=1, positiveStrand=this.positiveStrand, cigar=Cigar.empty, range=range)
  }
}

object AlignmentBlock extends LazyLogging {

  private def fail(message: String) = throw new IllegalStateException(message)

  /** Builds an alignment block from a [[SamRecord]]
   *
   * @param rec the mapped record
   */
  def apply(rec: SamRecord): AlignmentBlock = {
    require(rec.mapped)
    val range = GenomicRange(refIndex = rec.refIndex, start = rec.start, end = rec.end)
    val cigar = rec.cigar
    val leadingClipping = cigar.iterator
      .takeWhile(_.operator.isClipping) // take leading clipping
      .map(_.length).sum
    val middle = cigar.iterator
      .dropWhile(_.operator.isClipping) // skip leading clipping
      .takeWhile(!_.operator.isClipping) // take until we hit any clipping at the end
      .filter(_.operator.consumesReadBases()) // only operators that consume read bases
      .map(_.length).sum
    val trailingClipping = cigar.iterator
      .dropWhile(_.operator.isClipping) // skip leading clipping
      .dropWhile(!_.operator.isClipping) // skip until we hit any clipping at the end
      .map(_.length).sum
    val (start, end) = {
      if (rec.positiveStrand) (leadingClipping + 1, leadingClipping + middle)
      else (trailingClipping + 1, trailingClipping + middle)
    }

    AlignmentBlock(origin = BlockOrigin(rec), start = start, end = end, positiveStrand = rec.positiveStrand, cigar = rec.cigar, range = range)
  }

  /** Builds [[AlignmentBlock]]s for the given alignments for a given read, one per record. The blocks returned
   * are ordered by the first aligned base in the read in sequencing order.
   *
   * See [[blocksFrom()]] for a detailed explanation.
   *
   * @param primary             the primary alignment
   * @param supplementals       the supplemental alignments for the read
   * @param minUniqueBasesToAdd the minimum number of unique read bases the block must cover when adding to the current
   *                            set of blocks.
   */
  def blocksFrom(primary: SamRecord,
                 supplementals: Iterator[SamRecord],
                 minUniqueBasesToAdd: Int): IndexedSeq[AlignmentBlock] = {
    require(supplementals.nonEmpty)
    val primaryBlock = AlignmentBlock(primary)
    val supplBlocks  = supplementals.map(AlignmentBlock(_))
    this.blocksFrom(primary=primaryBlock, supplementals=supplBlocks, readLength=primary.length, minUniqueBasesToAdd=minUniqueBasesToAdd)
  }


  /** Builds [[AlignmentBlock]]s for the given alignments for a given read, one per record. The blocks returned
   * are ordered by the first aligned base in the read in sequencing order.
   *
   * The primary alignment is always added first.  Next, the supplementary alignment blocks are iterated in order from
   * offset from the first base in sequencing order.  For each supplementary alignment block, it is added to the list to
   * return only if the block covers at least a minimum # of currently uncovered query bases.
   *
   * @param primary             the primary alignment
   * @param supplementals       the supplemental alignments for the read
   * @param readLength          the length of the read in bases
   * @param minUniqueBasesToAdd the minimum number of unique read bases the block must cover when adding to the current
   *                            set of blocks.
   */
  def blocksFrom(primary: AlignmentBlock,
                 supplementals: Iterator[AlignmentBlock],
                 readLength: Int,
                 minUniqueBasesToAdd: Int): IndexedSeq[AlignmentBlock] = {
    require(supplementals.nonEmpty)
    val supplBlocks = supplementals.toIndexedSeq.sortBy(b => (b.start, b.end))
    val builder     = IndexedSeq.newBuilder[AlignmentBlock]

    require(supplementals.forall(_.origin == primary.origin))

    // how do we pick blocks to keep
    val mask = new mutable.BitSet(initSize=readLength)

    // add the primary one
    mask.addAll(primary.toQueryBitSet)
    builder += primary

    // figure out which supps to add
    supplBlocks.foreach { block =>
      // Add it to the tree if the block covers at least a minimum # of currently uncovered query bases.
      val blockMask = block.toQueryBitSet
      val uncovered = (blockMask &~ mask).size // remove the current bits from the block bits
      if (uncovered >= minUniqueBasesToAdd) {
        mask |= blockMask // union the bit masks
        builder += block
      }
    }

    builder.result().sortBy(b => (b.start, b.end))
  }


  /** Builds alignment blocks from a template.
   *
   * If no supplementary alignments exist, then an index block for R1 and R1 respectively is returned.  Otherwise,
   * supplementary alignments are filtered to have a minimum mapping quality, blocks built for reach, and then the
   * following iterative algorithm is performed for each read end:
   * 1. The primary alignment for the given read end is added to the set to keep.
   * 2. Next, the alignment block for the supplementary alignment that maps the left-most base in sequencing order is
   * chosen.  If it maps at least the `minUniqueBasesToAdd` new bases versus those alignment blocks already picked,
   * it is added to the set.  This is performed until no more supplementary alignments remain.
   * 3. The blocks for R1 and R2 produced in step (2) are merged with the procedure described in [[mergeReadBlocks()]].
   *
   * @param template                       the template from which alignment blocks are to be produced
   * @param minSupplementaryMappingQuality the minimum mapping quality to keep a supplementary alignment
   * @param minUniqueBasesToAdd            the minimum # of new bases that a supplementary alignment must map to keep it in the
   *                                       iterative procedure described above.
   */
  def blocksFrom(template: Template,
                 minSupplementaryMappingQuality: Int,
                 minUniqueBasesToAdd: Int): IndexedSeq[AlignmentBlock] = {
    val r1 = template.r1.getOrElse(fail(s"No R1 for template: ${template.name}"))
    val r2 = template.r2.getOrElse(fail(s"No R2 for template: ${template.name}"))
    val r1Supps = template.r1Supplementals.iterator.filter(_.mapq >= minSupplementaryMappingQuality)
    val r2Supps = template.r2Supplementals.iterator.filter(_.mapq >= minSupplementaryMappingQuality)

    if (r1Supps.isEmpty && r2Supps.isEmpty) {
      IndexedSeq(AlignmentBlock(r1), AlignmentBlock(r2))
    } else {
      val r1Blocks = if (r1Supps.isEmpty) IndexedSeq(AlignmentBlock(r1)) else {
        AlignmentBlock.blocksFrom(primary=r1, supplementals=r1Supps, minUniqueBasesToAdd=minUniqueBasesToAdd)
      }
      val r2Blocks = if (r2Supps.isEmpty) IndexedSeq(AlignmentBlock(r2)) else {
        AlignmentBlock.blocksFrom(primary=r2, supplementals=r2Supps, minUniqueBasesToAdd=minUniqueBasesToAdd)
          .reverse
          .map(b => b.copy(positiveStrand = !b.positiveStrand))
      }
      // reverse the R2 blocks as the first blocks in R2 are last blocks in the template when starting from the start
      // of R1, also need to switch strand as we expect FR pair
      mergeReadBlocks(r1Blocks=r1Blocks, r2Blocks=r2Blocks, numOverlappingBlocks=1)
    }
  }

  /**
   * Merges alignment blocks from read one and two respectively.
   *
   * It is assumed that the read pair is FR in orientation, and that the read two blocks are in reverse sequencing
   * order, so that the blocks that map the end of R1 are at the end of the provided sequence, while the blocks that map
   * at the end of R2 are at the start of the provided sequence.
   *
   * The last `numOverlappingBlocks` of R1 and the first `numOverlappingBlocks` of R2 are pairwise compared.  If all
   * pairs of blocks overlap in the genome and are on the same strand, they are merged, to create a chain of blocks.
   * If not all match, then this method is recursively called with `numOverlappingBlocks + 1`.  This repeats until
   * `numOverlappingBlocks` exceeds the number of R1 or R2 blocks.
   *
   * @param r1Blocks             the alignment blocks from read one sorted ascending the start of the read in sequencing order.
   * @param r2Blocks             the alignment blocks from read two sorted ascending the start of the read in sequencing order.
   * @param numOverlappingBlocks the # of overlapping blocks to examine
   * @return
   */
  @tailrec
  def mergeReadBlocks(r1Blocks: IndexedSeq[AlignmentBlock],
                      r2Blocks: IndexedSeq[AlignmentBlock],
                      numOverlappingBlocks: Int): IndexedSeq[AlignmentBlock] = {
    if (numOverlappingBlocks > r1Blocks.length || numOverlappingBlocks > r2Blocks.length) {
      r1Blocks ++ r2Blocks // no overlapping blocks
    } else {
      val r1SubBlocks = r1Blocks.takeRight(numOverlappingBlocks)
      val r2SubBlocks = r2Blocks.take(numOverlappingBlocks)
      val overlaps    = r1SubBlocks.zip(r2SubBlocks).forall { case (b1, b2) => b1.overlapsWithStrand(b2) }
      if (overlaps) {
        // merge and return
        val leading  = r1Blocks.take(r1Blocks.length - numOverlappingBlocks)
        val trailing = r2Blocks.takeRight(r2Blocks.length - numOverlappingBlocks)
        val middle   = r1SubBlocks.zip(r2SubBlocks).map { case (b1, b2) => b1.merge(b2) }
        leading ++ middle ++ trailing
      } else {
        // recurse
        mergeReadBlocks(r1Blocks=r1Blocks, r2Blocks=r2Blocks, numOverlappingBlocks=numOverlappingBlocks + 1)
      }
    }
  }
}

sealed trait BlockOrigin extends EnumEntry {
  /** Returns whether this origin and the other origin are R1 and R2, or R2 and R1, respectively.
   * This treats Both as being from either R1 and R2. */
  def isPairedWith(other: BlockOrigin): Boolean = this != other || (this == BlockOrigin.Both && other == BlockOrigin.Both)
}

/** Enumeration for where a given [[AlignmentBlock]] origintes. */
object BlockOrigin extends FgBioEnum[BlockOrigin] {
  def values: IndexedSeq[BlockOrigin] = findValues
  /** R1-only block */
  case object ReadOne extends BlockOrigin
  /** R2-only block */
  case object ReadTwo extends BlockOrigin
  /** Block merged from overlapping split read mapping in R1 and R2 */
  case object Both extends BlockOrigin
  /** Gets the origin from a SAM record */
  def apply(rec: SamRecord): BlockOrigin = {
    if (rec.firstOfPair) ReadOne else ReadTwo
  }
}