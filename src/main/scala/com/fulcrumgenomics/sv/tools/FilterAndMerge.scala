package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.NumericCounter
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sv.BreakpointPileup
import com.fulcrumgenomics.util.{Io, Metric}

import scala.collection.immutable.IndexedSeq

@clp(group=ClpGroups.BreakpointAndSv, description=
  """
    |Filters and merges SVPileup output.
  """)
class FilterAndMerge
( @arg(flag='i', doc="The input pileup file from SvPileup") input: FilePath,
  @arg(flag='o', doc="The output filtered and merged SvPileup file") output: FilePath,
  @arg(flag='m', doc="The minimum # of observations to examine an input site") minPre: Int = 1,
  @arg(flag='M', doc="The minimum # of observations to output a site") minPost: Int = 1,
  @arg(flag='s', doc="The maximum # bases between a breakend across adjacent sites") slop: Int = 0,
) extends SvTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    // Read in the pileups
    val iter   = Metric.iterator[BreakpointPileup](input)
      .filter(_.total >= minPre)
      .toIndexedSeq
      .sortBy(p => (p.left_contig, p.right_contig, p.left_strand, p.right_strand, p.left_pos, p.right_pos))
      .iterator
      .bufferBetter

    // Iterate
    val results = IndexedSeq.newBuilder[MergedPileup]
    while (iter.hasNext) {
      val pileups = nextGroup(iter=iter)
      val merged  = MergedPileup(pileups)
      if (merged.total_evidence >= minPost) {
        results += merged
      }
    }

    // Write
    Metric.write(output, results.result().sorted)
  }

  private def nextGroup(iter: BetterBufferedIterator[BreakpointPileup]): IndexedSeq[BreakpointPileup] = {
    val first   = iter.head
    val builder = IndexedSeq.newBuilder[BreakpointPileup]
    while (iter.hasNext && nearby(first, iter.head, slop)) {
      builder += iter.next()
    }

    builder.result()
  }

  /** Are two breakpoints compatible with each other and have positions that are +/- slop on each breakend. */
  private def nearby(a: BreakpointPileup, b: BreakpointPileup, slop: Int): Boolean = {
    a.left_contig == b.left_contig &&
      a.right_contig == b.right_contig &&
      a.left_strand == b.left_strand &&
      a.right_strand == b.right_strand &&
      Math.abs(a.left_pos - b.left_pos) <= slop &&
      Math.abs(a.right_pos - b.right_pos) <= slop
  }
}

case class MergedPileup(left_contig: String,
                        left_start: Int,
                        left_end: Int,
                        left_mean: Int,
                        left_strand: Char,
                        right_contig: String,
                        right_start: Int,
                        right_end: Int,
                        right_mean: Int,
                        right_strand: Char,
                        mean_count: Double,
                        stddev_count: Double,
                        median_count: Double,
                        pileups: Int,
                        split_reads: Int,
                        read_pairs: Int,
                        total_evidence: Int
                       ) extends Metric with Ordered[MergedPileup] {
  override def compare(that: MergedPileup): Int = {
    var result              = this.left_contig.compare(that.left_contig)
    if (result == 0) result = this.right_contig.compare(that.right_contig)
    if (result == 0) result = this.left_start.compare(that.left_start)
    if (result == 0) result = this.right_start.compare(that.right_start)
    if (result == 0) result = this.left_strand.compare(that.left_strand)
    if (result == 0) result = this.right_strand.compare(that.right_strand)
    if (result == 0) result = this.total_evidence.compare(that.total_evidence)
    result
  }
}

object MergedPileup {
  def apply(pileups: IndexedSeq[BreakpointPileup]): MergedPileup = {
    val splitReads = pileups.sumBy(_.split_reads)
    val readPairs  = pileups.sumBy(_.read_pairs)
    val counts     = NumericCounter(pileups.map(_.total).iterator)
    val total      = splitReads + readPairs
    val leftMean   = pileups.iterator.map(p => p.total * p.left_pos / total).sum
    val rightMean  = pileups.iterator.map(p => p.total * p.right_pos / total).sum
    val meanCount  = counts.mean()

    new MergedPileup(
      left_contig    = pileups.head.left_contig,
      left_start     = pileups.map(_.left_pos).min,
      left_end       = pileups.map(_.left_pos).max,
      left_mean      = Math.round(leftMean.toFloat),
      left_strand    = pileups.head.left_strand,
      right_contig   = pileups.head.right_contig,
      right_start    = pileups.map(_.right_pos).min,
      right_end      = pileups.map(_.right_pos).max,
      right_mean     = Math.round(rightMean.toFloat),
      right_strand   = pileups.head.right_strand,
      mean_count     = meanCount,
      stddev_count   = counts.stddev(meanCount),
      median_count   = counts.median(),
      pileups        = pileups.length,
      split_reads    = splitReads,
      read_pairs     = readPairs,
      total_evidence = splitReads + readPairs
    )
  }
}