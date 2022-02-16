package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{PathPrefix, PathToBam, SafelyClosable, yieldAndThen}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.EvidenceType._
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sv.{AlignedSegment, Breakpoint, BreakpointEvidence, EvidenceType}
import com.fulcrumgenomics.util.{Io, ProgressLogger}

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable


@clp(group=ClpGroups.All, description=
  """
    |Collates a pileup of putative structural variant supporting reads.
    |
    |Two output files will be created:
    |
    |1. `<output-prefix>.txt`: a tab-delimited file describing SV pileups, one line per breakpiont event.
    |2. `<output-prefix>.bam`: a SAM/BAM file containing reads that contain SV breakpoint evidence annotated with SAM
    |  tags.  The `ev` SAM tag lists the type of evidence found, while the `be` list the unique breakpoint identifier
    |  output in (1) above.
  """)
class SvPileup
(@arg(flag='i', doc="The input query sorted or grouped BAM") input: PathToBam,
 @arg(flag='o', doc="The output path prefix") output: PathPrefix,
 @arg(flag='d', doc="The maximum _inner_ distance for normal read pair") maxReadPairInnerDistance: Int = 1000,
 @arg(flag='D', doc="The maximum _inner_ distance between two segments of a split read mapping") maxAlignedSegmentInnerDistance: Int = 100,
 @arg(flag='q', doc="The minimum mapping quality for primary alignments") minPrimaryMappingQuality: Int = 30,
 @arg(flag='Q', doc="The minimum mapping quality for supplementary alignments") minSupplementaryMappingQuality: Int = 18,
 @arg(flag='b', doc="The minimum # of uncovered query bases needed to add a supplemental alignment") minUniqueBasesToAdd: Int = 20,
) extends SvTool {

  import SvPileup._

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val source     = SamSource(PathUtil.pathTo(output + ".bam"))
    val writer     = SamWriter(PathUtil.pathTo(output + ".txt"), header=source.header)
    val progress   = new ProgressLogger(logger, noun="templates")
    val tracker    = new BreakpointTracker()
    Bams.templateIterator(source)
      .tapEach(t => progress.record(t.r1.getOrElse(t.r2.get)))
      .flatMap(template => filterTemplate(template, minPrimaryMapq=minPrimaryMappingQuality, minSupplementaryMapq=minSupplementaryMappingQuality))
      .foreach { template =>
        // Find the breakpoints
        val evidences = findBreakpoints(
          template                       = template,
          maxWithinReadDistance          = maxAlignedSegmentInnerDistance,
          maxReadPairInnerDistance       = maxReadPairInnerDistance,
          minUniqueBasesToAdd            = minUniqueBasesToAdd,
        )

        // Update the tracker
        evidences.foreach { ev => tracker.count(ev.breakpoint, ev.evidence) }

        // Optionally write the reads to a BAM
        maybeWriteTemplate(template, evidences, tracker, writer)
      }

    progress.logLast()
    source.safelyClose()
    writer.close()

    // Write the results
    writeBreakpoints(tracker=tracker, dict=source.dict)
  }

  /** Annotates the reads for the given template and writes them to the writer if provided.
   * */
  private def maybeWriteTemplate(template: Template,
                                 evidences: IndexedSeq[BreakpointEvidence],
                                 tracker: BreakpointTracker,
                                 writer: SamWriter): Unit = {
    if (evidences.nonEmpty) {
      val bps = evidences.map(e => tracker(e.breakpoint).id).toSet.toSeq.sorted.mkString(",")
      val evs = evidences.map(_.evidence.snakeName).mkString(",")
      template.allReads.foreach { rec =>
        rec(SamBreakpointTag) = bps
        rec(SamEvidenceTag)   = evs
        writer += rec
      }
    }
  }

  /** Coalesce the breakpoint counts and write them. */
  private def writeBreakpoints(tracker: BreakpointTracker, dict: SequenceDictionary): Unit = {
    val writer = Io.toWriter(output)
    val fields = Seq("id", "left_contig", "left_pos", "left_strand", "right_contig", "right_pos", "right_strand",
      "split_reads", "read_pairs", "total")
    writer.write(fields.mkString("", "\t", "\n"))

    val breakpoints = tracker.breakpoints
      .toIndexedSeq
      .sortWith(Breakpoint.PairedOrdering.lt)

    breakpoints.foreach { bp =>
      val info         = tracker(bp)
      val id           = info.id
      val leftRefName  = dict(bp.leftRefIndex).name
      val rightRefName = dict(bp.rightRefIndex).name

      val values = Iterator(
        id, leftRefName, bp.leftPos, toStrand(bp.leftPositive), rightRefName, bp.rightPos, toStrand(bp.rightPositive),
        info.splitRead, info.readPair, info.total
      )
      writer.write(values.mkString("", "\t", "\n"))
    }

    writer.close()
  }

  /** Converts the boolean strand info to +/-. */
  private def toStrand(positive: Boolean): String = if (positive) "+" else "-"
}

object SvPileup extends LazyLogging {
  val SamEvidenceTag: String = "ev"
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
   * Performs filtering on a Template to remove primary records that are unmapped and then remove
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
   */
  def findBreakpoints(template: Template,
                      maxWithinReadDistance: Int,
                      maxReadPairInnerDistance: Int,
                      minUniqueBasesToAdd: Int,
                     ): IndexedSeq[BreakpointEvidence] = {
    val segments = AlignedSegment.segmentsFrom(template, minUniqueBasesToAdd=minUniqueBasesToAdd)

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
      val bp = Breakpoint(seg1, seg2)
      val ev = if (seg1.origin.isInterRead(seg2.origin)) EvidenceType.ReadPair else EvidenceType.SplitRead
      Some(BreakpointEvidence(bp, ev))
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
