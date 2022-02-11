package com.fulcrumgenomics.sv

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.util.LazyLogging
import enumeratum.EnumEntry

import scala.annotation.tailrec
import scala.collection.immutable.IndexedSeq
import scala.collection.{BitSet, mutable}

/** Represents an alignment segment of read to the reference. The alignment may not consume all the read bases.
 *
 * @param origin         from which read (or both) this segment comes from
 * @param readStart      the 1-based start in the read **in sequencing order**
 * @param readEnd        the 1-based inclusive end in the read **in sequencing order**
 * @param positiveStrand true if the read was mapped to the positive strand, false otherwise
 * @param cigar          the cigar **in sequencing order**
 * @param range          the genomic range to which this alignment maps
 */
case class AlignedSegment(origin: SegmentOrigin,
                          private val readStart: Int,
                          private val readEnd: Int,
                          positiveStrand: Boolean,
                          private val cigar: Cigar,
                          range: GenomicRange) {
  require(0 < readStart)
  require(readStart <= readEnd)

  /** Returns a bit set where each bit is a query base (in sequencing order), and bit is set if this alignment maps a base */
  def toQueryBitSet: BitSet = {
    scala.collection.BitSet(Range.inclusive(start=this.readStart, end=this.readEnd): _*)
  }

  /** Returns if this alignment overlaps the other alignment on the genome, and maps to the same strand. */
  def overlapsWithStrand(other: AlignedSegment): Boolean = {
    this.range.overlaps(other.range) && this.positiveStrand == other.positiveStrand
  }

  /** Merges the other segment with itself, returning a new alignment segment.  The two segments must overlap on the genome
   * and map to teh same strand.  If merging segments from different reads (i.e. R1 and R2), the segment origin will be set
   * to both.  The resulting cigar will be empty, and the start and end offset in the read will both be 1.
   * */
  def merge(other: AlignedSegment): AlignedSegment = {
    require(this.overlapsWithStrand(other))
    val range = this.range.union(other.range)
    val origin = if (this.origin == other.origin) this.origin else SegmentOrigin.Both
    AlignedSegment(origin=origin, readStart=1, readEnd=1, positiveStrand=this.positiveStrand, cigar=Cigar.empty, range=range)
  }
}

object AlignedSegment extends LazyLogging {

  private def fail(message: String) = throw new IllegalStateException(message)

  /** Builds an alignment segment from a [[SamRecord]]
   *
   * @param rec the mapped record
   */
  def apply(rec: SamRecord): AlignedSegment = {
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

    AlignedSegment(origin = SegmentOrigin(rec), readStart = start, readEnd = end, positiveStrand = rec.positiveStrand, cigar = rec.cigar, range = range)
  }

  /** Builds [[AlignedSegment]]s for the given alignments for a given read, one per record. The segments returned
   * are ordered by the first aligned base in the read in sequencing order.
   *
   * See [[segmentsFrom()]] for a detailed explanation.
   *
   * @param primary             the primary alignment
   * @param supplementals       the supplemental alignments for the read
   * @param minUniqueBasesToAdd the minimum number of unique read bases the segment must cover when adding to the current
   *                            set of segments.
   */
  def segmentsFrom(primary: SamRecord,
                   supplementals: Iterator[SamRecord],
                   minUniqueBasesToAdd: Int): IndexedSeq[AlignedSegment] = {
    require(supplementals.nonEmpty)
    val primarySegment = AlignedSegment(primary)
    val supplSegments  = supplementals.map(AlignedSegment(_))
    this.segmentsFrom(primary=primarySegment, supplementals=supplSegments, readLength=primary.length, minUniqueBasesToAdd=minUniqueBasesToAdd)
  }


  /** Builds [[AlignedSegment]]s for the given alignments for a given read, one per record. The segments returned
   * are ordered by the first aligned base in the read in sequencing order.
   *
   * The primary alignment is always added first.  Next, the supplementary alignment segments are iterated in order from
   * offset from the first base in sequencing order.  For each supplementary alignment segment, it is added to the list to
   * return only if the segment covers at least a minimum # of currently uncovered query bases.
   *
   * @param primary             the primary alignment
   * @param supplementals       the supplemental alignments for the read
   * @param readLength          the length of the read in bases
   * @param minUniqueBasesToAdd the minimum number of unique read bases the segment must cover when adding to the current
   *                            set of segments.
   */
  def segmentsFrom(primary: AlignedSegment,
                   supplementals: Iterator[AlignedSegment],
                   readLength: Int,
                   minUniqueBasesToAdd: Int): IndexedSeq[AlignedSegment] = {
    require(supplementals.nonEmpty)
    val supplSegments = supplementals.toIndexedSeq.sortBy(b => (b.readStart, b.readEnd))
    val builder     = IndexedSeq.newBuilder[AlignedSegment]

    require(supplementals.forall(_.origin == primary.origin))

    // how do we pick segments to keep
    val mask = new mutable.BitSet(initSize=readLength)

    // add the primary one
    mask.addAll(primary.toQueryBitSet)
    builder += primary

    // figure out which supps to add
    supplSegments.foreach { segment =>
      // Add it to the tree if the segment covers at least a minimum # of currently uncovered query bases.
      val segmentMask = segment.toQueryBitSet
      val uncovered = (segmentMask &~ mask).size // remove the current bits from the segment bits
      if (uncovered >= minUniqueBasesToAdd) {
        mask |= segmentMask // union the bit masks
        builder += segment
      }
    }

    builder.result().sortBy(b => (b.readStart, b.readEnd))
  }


  /** Builds alignment segments from a template.
   *
   * If no supplementary alignments exist, then an index segment for R1 and R1 respectively is returned.  Otherwise,
   * supplementary alignments are filtered to have a minimum mapping quality, segments built for reach, and then the
   * following iterative algorithm is performed for each read end:
   * 1. The primary alignment for the given read end is added to the set to keep.
   * 2. Next, the alignment segment for the supplementary alignment that maps the left-most base in sequencing order is
   * chosen.  If it maps at least the `minUniqueBasesToAdd` new bases versus those alignment segments already picked,
   * it is added to the set.  This is performed until no more supplementary alignments remain.
   * 3. The segments for R1 and R2 produced in step (2) are merged with the procedure described in [[mergeAlignedSegments()]].
   *
   * @param template                       the template from which alignment segments are to be produced
   * @param minSupplementaryMappingQuality the minimum mapping quality to keep a supplementary alignment
   * @param minUniqueBasesToAdd            the minimum # of new bases that a supplementary alignment must map to keep it in the
   *                                       iterative procedure described above.
   */
  def segmentsFrom(template: Template,
                   minSupplementaryMappingQuality: Int,
                   minUniqueBasesToAdd: Int): IndexedSeq[AlignedSegment] = {
    val r1 = template.r1.getOrElse(fail(s"No R1 for template: ${template.name}"))
    val r2 = template.r2.getOrElse(fail(s"No R2 for template: ${template.name}"))
    val r1Supps = template.r1Supplementals.iterator.filter(_.mapq >= minSupplementaryMappingQuality)
    val r2Supps = template.r2Supplementals.iterator.filter(_.mapq >= minSupplementaryMappingQuality)

    if (r1Supps.isEmpty && r2Supps.isEmpty) {
      IndexedSeq(AlignedSegment(r1), AlignedSegment(r2))
    } else {
      val r1Segments = if (r1Supps.isEmpty) IndexedSeq(AlignedSegment(r1)) else {
        AlignedSegment.segmentsFrom(primary=r1, supplementals=r1Supps, minUniqueBasesToAdd=minUniqueBasesToAdd)
      }
      val r2Segments = if (r2Supps.isEmpty) IndexedSeq(AlignedSegment(r2)) else {
        AlignedSegment.segmentsFrom(primary=r2, supplementals=r2Supps, minUniqueBasesToAdd=minUniqueBasesToAdd)
          .reverse
          .map(b => b.copy(positiveStrand = !b.positiveStrand))
      }
      // reverse the R2 segments as the first segments in R2 are last segments in the template when starting from the start
      // of R1, also need to switch strand as we expect FR pair
      mergeAlignedSegments(r1Segments=r1Segments, r2Segments=r2Segments, numOverlappingSegments=1)
    }
  }

  /**
   * Merges alignment segments from read one and two respectively.
   *
   * It is assumed that the read pair is FR in orientation, and that the read two segments are in reverse sequencing
   * order, so that the segments that map the end of R1 are at the end of the provided sequence, while the segments that map
   * at the end of R2 are at the start of the provided sequence.
   *
   * The last `numOverlappingSegments` of R1 and the first `numOverlappingSegments` of R2 are pairwise compared.  If all
   * pairs of segments overlap in the genome and are on the same strand, they are merged, to create a chain of segments.
   * If not all match, then this method is recursively called with `numOverlappingSegments + 1`.  This repeats until
   * `numOverlappingSegments` exceeds the number of R1 or R2 segments.
   *
   * @param r1Segments             the alignment segments from read one sorted ascending the start of the read in sequencing order.
   * @param r2Segments             the alignment segments from read two sorted ascending the start of the read in sequencing order.
   * @param numOverlappingSegments the # of overlapping segments to examine
   * @return
   */
  @tailrec
  def mergeAlignedSegments(r1Segments: IndexedSeq[AlignedSegment],
                           r2Segments: IndexedSeq[AlignedSegment],
                           numOverlappingSegments: Int): IndexedSeq[AlignedSegment] = {
    if (numOverlappingSegments > r1Segments.length || numOverlappingSegments > r2Segments.length) {
      r1Segments ++ r2Segments // no overlapping segments
    } else {
      val r1SubSegments = r1Segments.takeRight(numOverlappingSegments)
      val r2SubSegments = r2Segments.take(numOverlappingSegments)
      val overlaps    = r1SubSegments.zip(r2SubSegments).forall { case (seg1, seg2) => seg1.overlapsWithStrand(seg2) }
      if (overlaps) {
        // merge and return
        val leading  = r1Segments.take(r1Segments.length - numOverlappingSegments)
        val trailing = r2Segments.takeRight(r2Segments.length - numOverlappingSegments)
        val middle   = r1SubSegments.zip(r2SubSegments).map { case (seg1, seg2) => seg1.merge(seg2) }
        leading ++ middle ++ trailing
      } else {
        // recurse
        mergeAlignedSegments(r1Segments=r1Segments, r2Segments=r2Segments, numOverlappingSegments=numOverlappingSegments + 1)
      }
    }
  }
}

sealed trait SegmentOrigin extends EnumEntry {
  /** Returns whether this origin and the other origin are R1 and R2, or R2 and R1, respectively.
   * This treats Both as being from either R1 and R2. */
  def isPairedWith(other: SegmentOrigin): Boolean = this != other || (this == SegmentOrigin.Both && other == SegmentOrigin.Both)
}

/** Enumeration for where a given [[AlignedSegment]] origintes. */
object SegmentOrigin extends FgBioEnum[SegmentOrigin] {
  def values: IndexedSeq[SegmentOrigin] = findValues
  /** R1-only segment */
  case object ReadOne extends SegmentOrigin
  /** R2-only segment */
  case object ReadTwo extends SegmentOrigin
  /** segment merged from overlapping split read mapping in R1 and R2 */
  case object Both extends SegmentOrigin
  /** Gets the origin from a SAM record */
  def apply(rec: SamRecord): SegmentOrigin = {
    if (rec.firstOfPair) ReadOne else ReadTwo
  }
}