package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, PathPrefix, PathToBam, SafelyClosable, yieldAndThen, _}
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.EvidenceType._
import com.fulcrumgenomics.sv._
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.util.{Io, Metric, ProgressLogger}
import enumeratum.EnumEntry
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.util.OverlapDetector
import htsjdk.tribble.bed.BEDFeature

import scala.collection.immutable.IndexedSeq
import scala.collection.{immutable, mutable}

/** An enumeration over how to utilize the input target BED file if given. */
sealed trait TargetBedRequirement extends EnumEntry

/** An enumeration over how to utilize the input target BED file if given. */
object TargetBedRequirement extends FgBioEnum[TargetBedRequirement] {
  /** Annotate the breakpionts only. */
  case object AnnotateOnly extends TargetBedRequirement
  /** Keep only breakpoints where at least one side overlaps at least one the targets. */
  case object OverlapAny extends TargetBedRequirement
  /** Keep only breakpoints where both sides overlap at least one of the targets. */
  case object OverlapBoth extends TargetBedRequirement

  override val values: immutable.IndexedSeq[TargetBedRequirement] = findValues

  type Requirement = TargetBedRequirement
}


@clp(group=ClpGroups.All, description=
  """
    |Collates a pileup of putative structural variant supporting reads.
    |
    |## Outputs
    |
    |Two output files will be created:
    |
    |1. `<output-prefix>.txt`: a tab-delimited file describing SV pileups, one line per breakpoint event.  The returned
    |   breakpoint will be canonicalized such that the "left" side of the breakpoint will have the lower (or equal to)
    |   position on the genome vs. the "right"s side.
    |2. `<output-prefix>.bam`: a SAM/BAM file containing reads that contain SV breakpoint evidence annotated with SAM
    |  tag.
    |
    |The `be` SAM tag contains a comma-delimited list of breakpoints to which a given read belongs.  Each element is
    |a semi-colon delimited, with four fields:
    |
    |1. The unique breakpoint identifier (same identifier found in the tab-delimited output).
    |2. Either "left" or "right, corresponding to if the read shows evidence of the genomic left or right side of the
    |   breakpoint as found in the breakpoint file (i.e. `left_pos` or `right_pos`).
    |3. Either "from" or "into", such that when traversing the breakpoint would read through "from" and then into
    |   "into" in the sequencing order of the read pair.  For a split-read alignment, the "from" contains the aligned
    |   portion of the read that comes from earlier in the read in sequencing order.  For an alignment of a read-pair
    |   spanning the breakpoint, then "from" should be read-one of the pair and "into" should be read-two of the pair.
    |4. The type of breakpoint evidence: either "split_read" for observations of an aligned segment of a single read
    |   with split alignments, or "read_pair" for observations _between_ reads in a read pair.
    |
    |## Example output
    |
    |The following shows two breakpoints:
    |
    |```
    |id left_contig left_pos left_strand right_contig right_pos right_strand split_reads read_pairs total
    | 1        chr1      100           +         chr2       200            -           1          0     1
    | 2        chr2      150           -         chr3       500            +           1          0     1
    |```
    |
    |Consider a single fragment read that maps across both the above two breakpoints, so has three split-read
    |alignments.  The first alignment maps on the left side of breakpoint #1, the second alignment maps to both the
    |right side of breakpoint #1 and the left-side of breakpoint #2, and the third alignment maps to the right side of
    |breakpoint #2. The SAM records would be as follows:
    |
    |```
    |r1    0 chr1  50 60   50M100S ... be:Z:1;left;from;split_read
    |r1 2064 chr2 150 60 50S50M50S ... be:Z:1;right;into;split_read,2;left;from;split_read
    |r1 2048 chr3 500 60   100S50M ... be:Z:2;right;into;split_read
    |```
    |
    |## Algorithm Overview
    |
    |Putative breakpoints are identified by examining the alignments for each template. The alignments are transformed
    |into aligned segments in the order they were sequenced.  Each aligned segment represents the full genomic span of
    |the mapped bases.  This is performed first for the primary alignments.  Next, supplementary alignments are added
    |only if they map read bases that have not been previously covered by other alignments (see
    |`--min-unique-bases-to-add`).  This is iteratively performed until supplementary alignments have been exhausted.
    |
    |Next, aligned segments that have overlapping genomic mapped bases are merged into a single aligned
    |segment.  In this case, the two or more read mappings merged are associated with either the left side or right
    |side of that aligned segment, controlled by examining how close to the end of the new aligned segment the given
    |read mapping occurs (see `--slop` option).  This used to identify which reads traverse "from" and "into" the
    |breakpoint as described above.
    |
    |Finally, pairs of adjacent aligned segments are examined for evidence of a breakpoint, such that genomic distance
    |between them beyond either `--max-read-pair-inner-distance` for aligned segments from different read pairs, or
    |`--max-aligned-segment-inner-distance` for aligned segments from the same read in a pair (i.e. split-read mapping).
    |Split read evidence will be returned in favor of across-read-pair evidence when both are present.
  """)
class SvPileup
(@arg(flag='i', doc="The input query sorted or grouped BAM") input: PathToBam,
 @arg(flag='o', doc="The output path prefix") output: PathPrefix,
 @arg(flag='d', doc="The maximum _inner_ distance for normal read pair") maxReadPairInnerDistance: Int = 1000,
 @arg(flag='D', doc="The maximum _inner_ distance between two segments of a split read mapping") maxAlignedSegmentInnerDistance: Int = 100,
 @arg(flag='q', doc="The minimum mapping quality for primary alignments") minPrimaryMappingQuality: Int = 30,
 @arg(flag='Q', doc="The minimum mapping quality for supplementary alignments") minSupplementaryMappingQuality: Int = 18,
 @arg(flag='b', doc="The minimum # of uncovered query bases needed to add a supplemental alignment") minUniqueBasesToAdd: Int = 20,
 @arg(flag='s', doc="""
   |The number of bases of slop to allow when determining which records to track for the left or right
   |side of an aligned segment when merging segments.""") slop: Int = 5,
 @arg(flag='t', doc="Optional bed file of target regions") targetsBed: Option[FilePath] = None,
 @arg(flag='T', doc="Requirement on if each side of the breakpoint must overlap a target.  Will always annotate each side of the breakpoint.")
 targetsBedRequirement: TargetBedRequirement.Requirement = TargetBedRequirement.AnnotateOnly
) extends SvTool {

  import SvPileup._

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val source    = SamSource(input)
    val outHeader = {
      val h = source.header.clone()
      h.setSortOrder(SortOrder.unsorted)
      h.setGroupOrder(GroupOrder.query)
      h
    }
    val bamOut    = SamWriter(PathUtil.pathTo(s"$output.bam"), header=outHeader)
    val progress  = new ProgressLogger(logger, noun="templates")
    val tracker   = new BreakpointTracker()
    val targets   = targetsBed.map(Intervals.overlapDetectorFrom)

    Bams.templateIterator(source)
      .tapEach(t => progress.record(t.allReads.next()))
      .flatMap(template => filterTemplate(template, minPrimaryMapq=minPrimaryMappingQuality, minSupplementaryMapq=minSupplementaryMappingQuality))
      .foreach { template =>
        // Find the breakpoints
        val evidences = findBreakpoints(
          template                 = template,
          maxWithinReadDistance    = maxAlignedSegmentInnerDistance,
          maxReadPairInnerDistance = maxReadPairInnerDistance,
          minUniqueBasesToAdd      = minUniqueBasesToAdd,
          slop                     = slop
        )

        val filteredEvidences = targets match {
          case None           => evidences
          case Some(detector) =>
            evidences.filter { ev =>
              val leftBreakpoint = ev.breakpoint.leftInterval(source.dict)
              val rightBreakpoint = ev.breakpoint.rightInterval(source.dict)
              targetsBedRequirement match {
                case TargetBedRequirement.AnnotateOnly => true
                case TargetBedRequirement.OverlapAny   => detector.overlapsAny(leftBreakpoint) || detector.overlapsAny(rightBreakpoint)
                case TargetBedRequirement.OverlapBoth  => detector.overlapsAny(leftBreakpoint) && detector.overlapsAny(rightBreakpoint)
              }
            }
        }

        // Update the tracker
        filteredEvidences.foreach { ev => tracker.count(ev.breakpoint, ev.evidence) }

        // Optionally write the reads to a BAM
        maybeWriteTemplate(template, filteredEvidences, tracker, bamOut)
      }

    progress.logLast()
    source.safelyClose()
    bamOut.close()

    // Write the results
    writeBreakpoints(path=PathUtil.pathTo(s"$output.txt"), tracker=tracker, dict=source.dict, targets=targets)
  }

  /** Annotates the reads for the given template and writes them to the writer if provided.
   * */
  private def maybeWriteTemplate(template: Template,
                                 evidences: IndexedSeq[BreakpointEvidence],
                                 tracker: BreakpointTracker,
                                 writer: SamWriter): Unit = {
    if (evidences.nonEmpty) {
      template.allReads.foreach { rec =>
        val builder = IndexedSeq.newBuilder[String]
        evidences.foreach { ev =>
          val id = tracker(ev.breakpoint).id
          if (ev.from.contains(rec)) {
            val leftOrRight = if (ev.fromIsLeft) "left" else "right"
            builder.addOne(f"$id;$leftOrRight;from;${ev.evidence.snakeName}")
          }
          if (ev.into.contains(rec)) {
            val leftOrRight = if (ev.fromIsLeft) "right" else "left"
            builder.addOne(f"$id;$leftOrRight;into;${ev.evidence.snakeName}")
          }
        }
        val values = builder.result()
        if (values.nonEmpty) rec(SamBreakpointTag) = values.mkString(",")
        writer += rec
      }
    }
  }

  /** Write the breakpoints to a file. */
  private def writeBreakpoints(path: FilePath,
                               tracker: BreakpointTracker,
                               dict: SequenceDictionary,
                               targets: Option[OverlapDetector[BEDFeature]]): Unit = {
    val writer      = Metric.writer[BreakpointPileup](path)
    val breakpoints = tracker.breakpoints.toIndexedSeq.sortWith(Breakpoint.PairedOrdering.lt)

    breakpoints.foreach { bp =>
      val leftTargets  = targets.map(_.getOverlaps(bp.leftInterval(dict)).map(_.getName).toSeq.sorted.distinct.mkString(",")).flatMap(s => if (s.isEmpty) None else Some(s))
      val rightTargets = targets.map(_.getOverlaps(bp.rightInterval(dict)).map(_.getName).toSeq.sorted.distinct.mkString(",")).flatMap(s => if (s.isEmpty) None else Some(s))
      val info         = tracker(bp)
      val metric       = BreakpointPileup(
        id            = info.id.toString,
        left_contig   = dict(bp.leftRefIndex).name,
        left_pos      = bp.leftPos,
        left_strand   = if (bp.leftPositive) '+' else '-',
        right_contig  = dict(bp.rightRefIndex).name,
        right_pos     = bp.rightPos,
        right_strand  = if (bp.rightPositive) '+' else '-',
        split_reads   = info.splitRead,
        read_pairs    = info.readPair,
        total         = info.total,
        left_targets  = leftTargets,
        right_targets = rightTargets
      )

      writer.write(metric)
    }

    writer.close()
  }
}

object SvPileup extends LazyLogging {
  val SamBreakpointTag: String = "be"

  type BreakpointId = Long

  /** Value used when no breakpoints are detected. */
  private val NoBreakpoints: IndexedSeq[BreakpointEvidence] = IndexedSeq.empty

  /** Tracks counts of split reads vs. read pairs supporting a Breakpoint. */
  private case class BreakpointInfo(id: BreakpointId, var splitRead: Int = 0 , var readPair: Int = 0) {
    def total: Int = splitRead + readPair
  }

  /** Class that tracks IDs and counts for Breakpoints during discovery. */
  private class BreakpointTracker extends Iterable[(Breakpoint, BreakpointInfo)] {
    private var _nextId: Long = 0L
    private val bpToInfo      = mutable.HashMap[Breakpoint, BreakpointInfo]()

    private def nextId: BreakpointId = yieldAndThen(this._nextId) { this._nextId += 1 }

    /** Adds a count of one to the evidence for the given breakpoint under the evidence type given.
     * Returns the numerical ID of the breakpoint.
     */
    def count(bp: Breakpoint, ev: EvidenceType): BreakpointId = {
      val info = this.bpToInfo.getOrElseUpdate(bp, BreakpointInfo(nextId))

      ev match {
        case SplitRead => info.splitRead += 1
        case ReadPair  => info.readPair  += 1
      }

      info.id
    }

    override def iterator: Iterator[(Breakpoint, BreakpointInfo)] = this.bpToInfo.iterator

    /** Returns an iterator over the set of observed breakpoints, ordering is not predictable. */
    def breakpoints: Iterator[Breakpoint] = bpToInfo.keysIterator

    /** Gets the evidence counts for a given breakpoint in the same order as EvidenceType.values. */
    def apply(bp: Breakpoint): BreakpointInfo = this.bpToInfo(bp)
  }

  /**
   * Performs filtering on a [[Template]] to remove primary records that are unmapped and then remove
   * supplementary records if there is no matching primary record retained.  If neither primary
   * record is retained, returns None, else returns Some(Template).
   */
  def filterTemplate(t: Template,
                     minPrimaryMapq: Int,
                     minSupplementaryMapq: Int): Option[Template] = {
    val r1PrimaryOk = t.r1.exists(r => r.mapped && r.mapq >= minPrimaryMapq)
    val r2PrimaryOk = t.r2.exists(r => r.mapped && r.mapq >= minPrimaryMapq)

    if (!r1PrimaryOk && !r2PrimaryOk) None else Some(
      Template(
        r1              = if (r1PrimaryOk) t.r1 else None,
        r2              = if (r2PrimaryOk) t.r2 else None,
        r1Supplementals = if (r1PrimaryOk) t.r1Supplementals.filter(_.mapq >= minSupplementaryMapq) else Nil,
        r2Supplementals = if (r2PrimaryOk) t.r2Supplementals.filter(_.mapq >= minSupplementaryMapq) else Nil,
      )
    )
  }


  /** Finds the breakpoints for the given template.
   *
   * @param template the template to examine
   * @param maxWithinReadDistance the maximum distance between two adjacent split read mappings before calling a breakpoint
   * @param maxReadPairInnerDistance the maximum inner distance between R1 and R2 before calling a breakpoint
   * @param minUniqueBasesToAdd the minimum newly covered bases to keep a supplementary alignment when iteratively
   *                            adding them.
   * @param slop                the number of bases of slop to allow when determining which records to track for the
   *                            left or right side of an aligned segment when merging segments
   */
  def findBreakpoints(template: Template,
                      maxWithinReadDistance: Int,
                      maxReadPairInnerDistance: Int,
                      minUniqueBasesToAdd: Int,
                      slop: Int = 0
                     ): IndexedSeq[BreakpointEvidence] = {
    val segments = AlignedSegment.segmentsFrom(template, minUniqueBasesToAdd=minUniqueBasesToAdd, slop=slop)

    segments.length match {
      case 0 | 1 =>
        NoBreakpoints
      case 2     =>
        // Special case for 2 since most templates will generate two segments and we'd like it to be efficient
        val bp = findBreakpoint(segments.head, segments.last, maxWithinReadDistance, maxReadPairInnerDistance)
        if (bp.isEmpty) NoBreakpoints else bp.toIndexedSeq
      case _     =>
        segments.iterator.sliding(2).flatMap { case Seq(seg1, seg2) =>
          findBreakpoint(seg1, seg2, maxWithinReadDistance, maxReadPairInnerDistance)
        }.toIndexedSeq
    }
  }

  /** Checks to see if there is a breakpoint between two segments and returns it, or None. */
  private def findBreakpoint(seg1: AlignedSegment,
                             seg2: AlignedSegment,
                             maxWithinReadDistance: Int,
                             maxReadPairInnerDistance: Int): Option[BreakpointEvidence] = {
    if (isInterContigBreakpoint(seg1, seg2) ||
        isIntraContigBreakpoint(seg1, seg2, maxWithinReadDistance, maxReadPairInnerDistance)
    ) {
      val ev = if (seg1.origin.isInterRead(seg2.origin)) EvidenceType.ReadPair else EvidenceType.SplitRead
      Some(BreakpointEvidence(from=seg1, into=seg2, evidence=ev))
    }
    else {
      None
    }
  }

  /** Determines if two segments are evidence of a breakpoint joining two different contigs.
   *
   * @param seg1 the first alignment segment
   * @param seg2 the second alignment segment
   */
  def isInterContigBreakpoint(seg1: AlignedSegment, seg2: AlignedSegment): Boolean = {
    val r1 = seg1.range
    val r2 = seg2.range
    r1.refIndex != r2.refIndex
  }

  /** Determines if the two segments are provide evidence of a breakpoint joining two different regions from
   * the same contig. Returns true if:
   *   - the two segments overlap (implying some kind of duplication) (note overlapping reads will get a merged seg)
   *   - the strand of the two segments differ (implying an inversion or other rearrangement)
   *   - the second segment is before the first segment on the genome
   *   - the distance between the two segments is larger than the maximum allowed (likely a deletion)
   *
   * @param seg1 the first alignment segment
   * @param seg2 the second alignment segment
   * @param maxWithinReadDistance the maximum distance between segments if they are from the same read
   * @param maxBetweenReadDistance the maximum distance between segments if they are from different reads
   */
  def isIntraContigBreakpoint(seg1: AlignedSegment,
                              seg2: AlignedSegment,
                              maxWithinReadDistance: Int,
                              maxBetweenReadDistance: Int): Boolean = {
    require(seg1.range.refIndex == seg2.range.refIndex)

    // The way aligned segments are generated for a template, if we have all the reads in the expected orientation
    // the segments should all come out on the same strand.  Therefore any difference in strand is odd.  In addition
    // any segment that "moves backwards" down the genome is odd, as genome position and read position should increase
    // together.
    if (seg1.positiveStrand != seg2.positiveStrand) true
    else if (seg1.positiveStrand && seg2.range.start < seg1.range.end) true
    else if (!seg1.positiveStrand && seg1.range.start < seg2.range.start) true
    else {
      val maxDistance = if (seg1.origin.isInterRead(seg2.origin)) maxBetweenReadDistance else maxWithinReadDistance

      val innerDistance = {
        if (seg1.range.start <= seg2.range.start) seg2.range.start - seg1.range.end
        else                                      seg1.range.start - seg2.range.end
      }

      innerDistance > maxDistance
    }
  }
}
