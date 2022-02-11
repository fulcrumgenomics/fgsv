package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{BetterBufferedIteratorScalaWrapper, FilePath, PathPrefix, PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sv.util.EvidenceType._
import com.fulcrumgenomics.sv.util.{AlignedSegment, EvidenceType, PutativeBreakpoint}
import com.fulcrumgenomics.sopt.{arg, clp}
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
    val source         = SamSource(PathUtil.pathTo(output + ".bam"))
    val writer         = SamWriter(PathUtil.pathTo(output + ".txt"), header=source.header)
    val progress       = new ProgressLogger(logger, noun="templates")
    val breakpoints    = new SimpleCounter[PutativeBreakpoint]()
    val breakpointToId = new BreakpointToIdMap()
    Bams.templateIterator(source)
      .map(t => {progress.record(t.r1.getOrElse(t.r2.get)); t})
      .filter(template => template.primaryReads.forall(primaryOk))
      .foreach { template =>
        // Find the breakpoints
        val templateBreakpoints = findBreakpoints(
          template                       = template,
          maxReadPairInnerDistance       = maxReadPairInnerDistance,
          maxAlignedSegmentInnerDistance = maxAlignedSegmentInnerDistance,
          minSupplementaryMappingQuality = minSupplementaryMappingQuality,
          minUniqueBasesToAdd            = minUniqueBasesToAdd,
          breakpointToId                 = breakpointToId
        )
        // Optionally write the reads to a BAM
        maybeWriteTemplate(
          template       = template,
          breakpoints    = templateBreakpoints,
          writer         = writer
        )
        // Count the breakpoints
        templateBreakpoints.foreach(breakpoints.count)
      }
    progress.logLast()
    source.safelyClose()
    writer.close()

    // Write the results
    writeBreakpoints(breakpoints=breakpoints, dict=source.dict)
  }

  /** Annotates the reads for the given template and writes them to the writer if provided.
   * */
  private def maybeWriteTemplate(template: Template,
                                 breakpoints: IndexedSeq[PutativeBreakpoint],
                                 writer: SamWriter): Unit = {
    if (breakpoints.nonEmpty) {
      val tagValue = breakpoints.map(_.id).toSet.toSeq.sorted.map(_.toString).mkString(",")
      val evidences = breakpoints.map(_.evidence.snakeName).mkString(",")
      template.allReads.foreach { rec =>
        rec(SamBreakpointTag) = tagValue
        rec(SamEvidenceTag)   = evidences
        writer += rec
      }
    }
  }

  /** Coalesce the breakpoint counts and write them. */
  private def writeBreakpoints(breakpoints: SimpleCounter[PutativeBreakpoint], dict: SequenceDictionary): Unit = {
    val writer = Io.toWriter(output)
    val fields = Seq("id", "left_contig", "left_pos", "right_contig", "right_pos", "same_strand", "total") ++
      EvidenceType.values.map(_.snakeName)
    writer.write(fields.mkString("", "\t", "\n"))
    writer.write('\n')
    val breakpointsIter = breakpoints.toIndexedSeq
      .sortBy { case (b, _) => (b.leftRefIndex, b.rightRefIndex, b.leftPos, b.rightPos, b.sameStrand)}
      .iterator.bufferBetter
    while (breakpointsIter.hasNext) {
      val first: PutativeBreakpoint = breakpointsIter.head._1
      val coalesced = breakpointsIter.takeWhile(_._1.sameBreakEnds(first)).toIndexedSeq
      val countsMap: Map[EvidenceType, Long] = coalesced
        .groupBy(_._1.evidence)
        .map { case (svType: EvidenceType, data: Seq[(PutativeBreakpoint, Long)]) => (svType, data.map(_._2).sum) }
      val leftRefName  = dict(first.leftRefIndex).name
      val rightRefName = dict(first.rightRefIndex).name
      val counts       = EvidenceType.values.map(svType => countsMap.getOrElse(svType, 0L))
      val total        = counts.sum
      // FIXME: wouldn't it be nice if this wasn't a comma-delimited list?  That would take a bunch of work to use
      // .sameBreakends and store evidence counts when initially gathering the breakpoints.  Defer to later.
      val id = coalesced.iterator.map(_._1.id).toSet.toSeq.sorted.mkString(",")
      writer.write(f"$id\t$leftRefName\t${first.leftPos}\t$rightRefName\t${first.rightPos}\t${first.sameStrand}\t$total")
      counts.foreach(count => writer.write(f"\t$count"))
      writer.write('\n')
    }
    writer.close()
  }

  /** Returns true if the primary alignment is mapped and has sufficient mapping quality */
  private def primaryOk(primary: SamRecord): Boolean = {
    primary.mapped && primary.mapq >= minPrimaryMappingQuality
  }
}

object SvPileup extends LazyLogging {

  val SamEvidenceTag: String = "ev"
  val SamBreakpointTag: String = "be"

  /** A little class to map a breakpoint to id, ignoring the existing id of the breakpoint. */
  private class BreakpointToIdMap {
    private var _nextId: Long = 0L
    private val breakpointToId = mutable.HashMap[PutativeBreakpoint, Long]()

    /** Returns the identifier of the breakpoint, if it exists, or the next identifier that should be used. */
    def apply(breakpoint: PutativeBreakpoint): Long = {
      // NB: use id=0 and evidence=any so we can find it in the map regardless of existing id and evidence type
      this.breakpointToId.getOrElseUpdate(breakpoint.copy(id=0, evidence=EvidenceType.ReadPairTandem), {
        val id = this._nextId
        this._nextId += 1
        id
      })
    }

    def nextId: Long = this._nextId
  }

  /** Finds the breakpoints for the given template.
   *
   * @param template the template to examine
   * @param maxReadPairInnerDistance the maximum inner distance between R1 and R2 before calling a breakpoint
   * @param maxAlignedSegmentInnerDistance the maximum distance between two adjacent split read mappings before calling
   *                                       a breakpoint
   * @param minSupplementaryMappingQuality the minimum mapping quailty to keep a supplementary alignment
   * @param minUniqueBasesToAdd the minimum newly covered bases to keep a supplementary alignment when iteratively
   *                            adding them.
   */
  def findBreakpoints(template: Template,
                      maxReadPairInnerDistance: Int,
                      maxAlignedSegmentInnerDistance: Int,
                      minSupplementaryMappingQuality: Int,
                      minUniqueBasesToAdd: Int,
                      breakpointToId: BreakpointToIdMap,
                     ): IndexedSeq[PutativeBreakpoint] = {
    val segments = AlignedSegment.segmentsFrom(
      template                       = template,
      minSupplementaryMappingQuality = minSupplementaryMappingQuality,
      minUniqueBasesToAdd            = minUniqueBasesToAdd
    )

    if (segments.length == 1) IndexedSeq.empty else {
      segments.iterator.sliding(2).flatMap { case Seq(seg1, seg2) =>
        var result: Option[PutativeBreakpoint] = determineInterContigBreakpiont(seg1=seg1, seg2=seg2)
        if (result.isEmpty) {
          val maxInnerDistance = if (seg1.origin != seg2.origin) maxAlignedSegmentInnerDistance else maxReadPairInnerDistance
          result = determineIntraContigBreakpiont(seg1=seg1, seg2=seg2, maxInnerDistance=maxInnerDistance)
        }
        if (result.isEmpty) result = determineOddPairOrientation(seg1=seg1, seg2=seg2)
        result
      }
      .map { breakpoint => breakpoint.copy(id = breakpointToId(breakpoint)) } // update the id!
      .toIndexedSeq
    }
  }

  /** Determines if two segments are evidence of an SV joining two different contigs.
   *
   * @param seg1 the first alignment segment
   * @param seg2 the second alignment segment
   */
  def determineInterContigBreakpiont(seg1: AlignedSegment, seg2: AlignedSegment): Option[PutativeBreakpoint] = {
    val r1 = seg1.range
    val r2 = seg2.range
    if (r1.refIndex == r2.refIndex) None
    else {
      val evidence = if (seg1.origin.isPairedWith(seg2.origin)) ReadPairInterContig else SplitReadInterContig
      Some(PutativeBreakpoint(seg1=seg1, seg2=seg2, evidence=evidence))
    }
  }

  /** Determines if the primary mappings are from a read pair are evidence of an SV joining two different regions from
   * the same contig.
   *
   * @param seg1 the first alignment segment
   * @param seg2 the second alignment segment
   * @param maxInnerDistance the maximum "inner distance" allowed between the two segments, otherwise return a putative
   *                         breakpoint
   */
  def determineIntraContigBreakpiont(seg1: AlignedSegment, seg2: AlignedSegment, maxInnerDistance: Int): Option[PutativeBreakpoint] = {
    require(seg1.range.refIndex == seg2.range.refIndex)
    val innerDistance = {
      if (seg1.range.start <= seg2.range.start) seg2.range.start - seg1.range.end
      else seg1.range.start - seg2.range.end
    }
    if (innerDistance <= maxInnerDistance) None else {
      val evidence = if (seg1.origin.isPairedWith(seg2.origin)) ReadPairIntraContig else SplitReadIntraContig
      Some(PutativeBreakpoint(seg1=seg1, seg2=seg2, evidence=evidence))
    }
  }

  /** Determines if the primary mappings are from a read pair are evidence of an inversion or other re-arrangement from
   * an unexpected (i.e. not FR) read pair orientation.  This only makes sense if
   * the segments are from the R1 and R2 primary alignments.
   *
   * @param seg1 the first alignment segment
   * @param seg2 the second alignment segment
   */
  def determineOddPairOrientation(seg1: AlignedSegment, seg2: AlignedSegment): Option[PutativeBreakpoint] = {
    require(seg1.range.refIndex == seg2.range.refIndex)

    if (seg1.origin.isPairedWith(seg2.origin)) { // Treat as R1->R2 or R2->R1
      if (seg1.positiveStrand == seg2.positiveStrand) Some(PutativeBreakpoint(seg1=seg1, seg2=seg2, evidence=ReadPairTandem))
      else {
        // Follows the implementation in htsjdk.samtool.SamPairUtil.getPairOrientation
        val (positiveStrandFivePrime, negativeStrandFivePrime) = {
          if (seg1.positiveStrand) (seg1.range.start, seg2.range.end) else (seg2.range.start, seg1.range.end)
        }
        if (positiveStrandFivePrime < negativeStrandFivePrime) None // FR is ok
        else Some(PutativeBreakpoint(seg1=seg1, seg2=seg2, evidence=ReadPairReverseForward)) // RF
      }
    }
    else {
      if (seg1.positiveStrand == seg2.positiveStrand) None
      else Some(PutativeBreakpoint(seg1=seg1, seg2=seg2, evidence=SplitReadOppositeStrand))
    }
  }
}



