package com.fulcrumgenomics.fgsv.tools

import com.fulcrumgenomics.FgBioDef.{BetterBufferedIteratorScalaWrapper, FilePath, PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.bam.{Bams, Template}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.fgsv.cmdline.{ClpGroups, SVTool}
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.fgsv.util.EvidenceType._
import com.fulcrumgenomics.fgsv.util.{AlignmentBlock, BlockOrigin, EvidenceType, PutativeBreakpoint}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable


@clp(group=ClpGroups.All, description=
  """
    |Collates a pileup of putative structural variant calls.
  """)
class SVPileup
( @arg(flag='i', doc="The input query sorted or grouped BAM") input: PathToBam,
  @arg(flag='o', doc="The output file") output: FilePath,
  @arg(flag='d', doc="The maximum _inner_ distance for normal read pair") maxReadPairInnerDistance: Int = 1000,
  @arg(flag='D', doc="The maximum _inner_ distance between two blocks of a split read mapping") maxAlignmentBlockInnerDistance: Int = 100,
  @arg(flag='q', doc="The minimum mapping quality for primary alignments") minPrimaryMappingQuality: Int = 30,
  @arg(flag='Q', doc="The minimum mapping quality for supplementary alignments") minSupplementaryMappingQuality: Int = 18,
  @arg(flag='b', doc="The minimum # of uncovered query bases needed to add a supplemental alignment") minUniqueBasesToAdd: Int = 20,
  @arg(doc="The path to the BAM file that support the breakpoints") bam: Option[PathToBam] = None,
  @arg(doc="The SAM tag to annotate records with `--bam`") tag: String = "be"
) extends SVTool {

  import SVPileup._

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  bam.foreach(Io.assertCanWriteFile(_))

  override def execute(): Unit = {
    val source         = SamSource(input)
    val writer         = bam.map(path => SamWriter(path, header=source.header, sort=Some(SamOrder.Queryname)))
    val progress       = new ProgressLogger(logger, noun="templates")
    val breakpoints    = new SimpleCounter[PutativeBreakpoint]()
    val breakpointToId = new BreakpointToIdMap()
    Bams.templateIterator(source)
      .map(t => {progress.record(t.r1.getOrElse(t.r2.get)); t})
      .filter(template => (template.r1 ++ template.r2).forall(primaryOk))
      .foreach { template =>
        // Find the breakpoints
        val templateBreakpoints = findBreakpoints(
          template                       = template,
          maxReadPairInnerDistance       = maxReadPairInnerDistance,
          maxAlignmentBlockInnerDistance = maxAlignmentBlockInnerDistance,
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
    writer.foreach(_.close())

    // Write the results
    writeBreakpoints(breakpoints=breakpoints, dict=source.dict)
  }

  /** Annotates the reads for the given template and writes them to the writer if provided.
   * */
  private def maybeWriteTemplate(template: Template,
                                 breakpoints: IndexedSeq[PutativeBreakpoint],
                                 writer: Option[SamWriter]): Unit = {
    if (breakpoints.nonEmpty) {
      writer.foreach { _writer =>
        val tagValue = breakpoints.map(_.id).toSet.toSeq.sorted.map(_.toString).mkString(",")
        val evidences = breakpoints.map(_.evidence.snakeName).mkString(",")
        template.allReads.foreach { rec =>
          rec(tag) = tagValue
          rec("ev") = evidences
          _writer += rec
        }
      }
    }
  }

  /** Coalesce the breakpoint counts and write them. */
  private def writeBreakpoints(breakpoints: SimpleCounter[PutativeBreakpoint], dict: SequenceDictionary): Unit = {
    val writer = Io.toWriter(output)
    writer.write("id\tleft_contig\tleft_pos\tright_contig\tright_pos\tsame_strand\ttotal")
    EvidenceType.values.foreach(svType => writer.write(s"\t${svType.snakeName}"))
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
      // .sameBreakends and store evidence couints when initially gathering the breakpoints.  Defer to later.
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

object SVPileup extends LazyLogging {

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
   * @param maxAlignmentBlockInnerDistance the maximum distance between two adjacent split read mappings before calling
   *                                       a breakpoint
   * @param minSupplementaryMappingQuality the minimum mapping quailty to keep a supplementary alignment
   * @param minUniqueBasesToAdd the minimum newly covered bases to keep a supplementary alignment when iteratively
   *                            adding them.
   */
  def findBreakpoints(template: Template,
                      maxReadPairInnerDistance: Int,
                      maxAlignmentBlockInnerDistance: Int,
                      minSupplementaryMappingQuality: Int,
                      minUniqueBasesToAdd: Int,
                      breakpointToId: BreakpointToIdMap,
                     ): IndexedSeq[PutativeBreakpoint] = {
    val blocks = AlignmentBlock.blocksFrom(
      template                       = template,
      minSupplementaryMappingQuality = minSupplementaryMappingQuality,
      minUniqueBasesToAdd            = minUniqueBasesToAdd
    )

    if (blocks.length == 1) IndexedSeq.empty else {
      blocks.iterator.sliding(2).flatMap { case Seq(b1, b2) =>
        var result: Option[PutativeBreakpoint] = determineInterContigChimera(b1=b1, b2=b2)
        if (result.isEmpty) {
          val maxInnerDistance = (b1.origin, b2.origin) match {
            case (BlockOrigin.ReadOne, BlockOrigin.ReadTwo) | (BlockOrigin.ReadTwo, BlockOrigin.ReadOne) => maxReadPairInnerDistance
            case _ => maxAlignmentBlockInnerDistance
          }
          result = determineIntraContigChimera(b1=b1, b2=b2, maxInnerDistance=maxInnerDistance)
        }
        if (result.isEmpty) result = determineOddPairOrientation(b1=b1, b2=b2)
        result
      }
      .map { breakpoint => breakpoint.copy(id = breakpointToId(breakpoint)) } // update the id!
      .toIndexedSeq
    }
  }

  /** Determines if the primary mappings are from a read pair are evidence of an SV joining two different contigs.
   *
   * @param b1 the first alignment block
   * @param b2 the second alignment block
   */
  def determineInterContigChimera(b1: AlignmentBlock, b2: AlignmentBlock): Option[PutativeBreakpoint] = {
    val r1 = b1.range
    val r2 = b2.range
    if (r1.refIndex == r2.refIndex) None
    else {
      val evidence = if (b1.origin.isPairedWith(b2.origin)) ReadPairInterContig else SplitReadInterContig
      Some(PutativeBreakpoint(b1=b1, b2=b2, evidence=evidence))
    }
  }

  /** Determines if the primary mappings are from a read pair are evidence of an SV joining two different regions from
   * the same contig.
   *
   * @param b1 the first alignment block
   * @param b2 the second alignment block
   * @param maxInnerDistance the maximum "inner distance" allowed between the two blocks, otherwise return a putative
   *                         breakpoint
   */
  def determineIntraContigChimera(b1: AlignmentBlock, b2: AlignmentBlock, maxInnerDistance: Int): Option[PutativeBreakpoint] = {
    require(b1.range.refIndex == b2.range.refIndex)
    val innerDistance = {
      if (b1.range.start <= b2.range.start) b2.range.start - b1.range.end
      else b1.range.start - b2.range.end
    }
    if (innerDistance <= maxInnerDistance) None else {
      val evidence = if (b1.origin.isPairedWith(b2.origin)) ReadPairIntraContig else SplitReadIntraContig
      Some(PutativeBreakpoint(b1=b1, b2=b2, evidence=evidence))
    }
  }

  /** Determines if the primary mappings are from a read pair are evidence of an inversion or other re-arrangement from
   * an unexpected (i.e. not FR) read pair orientation.  This only makes sense if
   * the blocks are from the R1 and R2 primary alignments.
   *
   * @param b1 the first alignment block
   * @param b2 the second alignment block
   */
  def determineOddPairOrientation(b1: AlignmentBlock, b2: AlignmentBlock): Option[PutativeBreakpoint] = {
    require(b1.range.refIndex == b2.range.refIndex)

    if (b1.origin.isPairedWith(b2.origin)) { // Treat as R1->R2 or R2->R1
      if (b1.positiveStrand == b2.positiveStrand) Some(PutativeBreakpoint(b1=b1, b2=b2, evidence=ReadPairTandem))
      else {
        // Follows the implementation in htsjdk.samtool.SamPairUtil.getPairOrientation
        val (positiveStrandFivePrime, negativeStrandFivePrime) = {
          if (b1.positiveStrand) (b1.start, b2.end) else (b2.start, b1.end)
        }
        if (positiveStrandFivePrime < negativeStrandFivePrime) None // FR is ok
        else Some(PutativeBreakpoint(b1=b1, b2=b2, evidence=ReadPairReverseForward)) // RF
      }
    }
    else {
      if (b1.positiveStrand == b2.positiveStrand) None
      else Some(PutativeBreakpoint(b1=b1, b2=b2, evidence=SplitReadOppositeStrand))
    }
  }
}



