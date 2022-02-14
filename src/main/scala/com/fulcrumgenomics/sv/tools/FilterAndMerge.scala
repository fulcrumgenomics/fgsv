package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{BetterBufferedIteratorScalaWrapper, FilePath, PathToSequenceDictionary}
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.{DelimitedDataParser, NumericCounter, Row, SimpleCounter}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.EvidenceType
import com.fulcrumgenomics.util.{Io, Metric}

import scala.annotation.tailrec
import scala.collection.immutable.IndexedSeq
import scala.collection.mutable

@clp(group=ClpGroups.All, description=
  """
    |Filters and merges SVPileup output.
  """)
class FilterAndMerge
( @arg(flag='i', doc="The input pileup file from SvPileup") input: FilePath,
  @arg(flag='o', doc="The output filtered and merged SvPileup file") output: FilePath,
  @arg(flag='d', doc="The path to the reference sequence dictionary.") dict: PathToSequenceDictionary,
  @arg(flag='m', doc="The minimum # of observations to examine an input site") minPre: Int = 1,
  @arg(flag='M', doc="The minimum # of observations to output a site") minPost: Int = 1,
  @arg(flag='s', doc="The maximum # bases between a breakend across adjacent sites") slop: Int = 0,
) extends SvTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val dict = SequenceDictionary(this.dict)

    // Read in the pileups
    val parser = DelimitedDataParser(path=input, delimiter='\t')
    val iter   = parser
      .map(row => Pileup(row))
      .filter(pileup => pileup.total >= minPre)
      .toIndexedSeq
      .groupBy(_.sameStrand) // this ensures we get pileups adjacent on the same strand
      .flatMap { case (_, pileups) =>
        pileups.sortBy(p => (p.leftContig, p.rightContig, p.leftPos, p.rightPos))
      }
      .iterator
      .bufferBetter

    // Iterate
    val results = IndexedSeq.newBuilder[MergedPileup]
    while (iter.hasNext) {
      val first   = iter.next()
      val builder = IndexedSeq.newBuilder[Pileup]
      builder += first
      val pileups = nextGroup(iter=iter, builder=builder, prev=first)
      val merged  = MergedPileup(pileups)
      if (merged.counter.total >= minPost) {
        results += merged
      }
    }

    // Write
    val writer = Io.toWriter(output)
    writer.write(MergedPileup.header.mkString("\t"))
    writer.write('\n')
    results.result()
      .sortBy { p =>
        (dict(p.left_contig).index, dict(p.right_contig).index, p.left_start, p.left_end, p.right_start, p.right_end)
      }
      .foreach { merged =>
        val values = MergedPileup.values(merged)
        writer.write(values.mkString("\t"))
        writer.write('\n')
      }
    writer.close()
  }

  @tailrec
  private def nextGroup(iter: BetterBufferedIterator[Pileup],
                        builder: mutable.Builder[Pileup, IndexedSeq[Pileup]],
                        prev: Pileup
                       ): IndexedSeq[Pileup] = {
    if (iter.isEmpty || !prev.nearby(other=iter.head, slop=slop)) {
      builder.result()
    } else {
      val next = iter.next()
      builder += next
      nextGroup(iter=iter, builder=builder, prev=next)
    }
  }
}

case class Pileup
(
  leftContig: String,
  leftPos: Int,
  rightContig: String,
  rightPos: Int,
  sameStrand: Boolean,
  counter: SimpleCounter[EvidenceType]
) {
  def total: Long = counter.total
  def nearby(other: Pileup, slop: Int): Boolean = {
    leftContig == other.leftContig &&
      rightContig == other.rightContig &&
      Math.abs(leftPos - other.leftPos) <= slop &&
      Math.abs(rightPos - other.rightPos) <= slop &&
      sameStrand == other.sameStrand
  }
}

object Pileup {
  def apply(row: Row): Pileup = {
    val counter = new SimpleCounter[EvidenceType]()
    EvidenceType.values.foreach { evidence =>
      counter.count(evidence, row[Int](evidence.snakeName))
    }
    new Pileup(
      leftContig  = row[String]("left_contig"),
      leftPos     = row[Int]("left_pos"),
      rightContig = row[String]("right_contig"),
      rightPos    = row[Int]("right_pos"),
      sameStrand  = row[Boolean]("same_strand"),
      counter     = counter
    )
  }
}

case class MergedPileup
(
  left_contig: String,
  left_start: Int,
  left_end: Int,
  left_mean: Int,
  right_contig: String,
  right_start: Int,
  right_end: Int,
  right_mean: Int,
  same_strand: Boolean,
  mean_count: Double,
  stddev_count: Double,
  median_count: Double,
  pileups: Int,
  counter: SimpleCounter[EvidenceType]
) extends Metric

object MergedPileup {
  def apply(pileups: IndexedSeq[Pileup]): MergedPileup = {
    val counter = new SimpleCounter[EvidenceType]()
    pileups.foreach { pileup =>
      counter += pileup.counter
    }

    val counts    = NumericCounter[Long](pileups.map(_.total).iterator)
    val total     = counts.totalMass.toDouble
    val leftMean  = pileups.iterator.map(p => p.total * p.leftPos / total).sum
    val rightMean = pileups.iterator.map(p => p.total * p.rightPos / total).sum
    val meanCount = counts.mean()
    new MergedPileup(
      left_contig  = pileups.head.leftContig,
      left_start   = pileups.map(_.leftPos).min,
      left_end     = pileups.map(_.leftPos).max,
      left_mean    = Math.round(leftMean).toInt,
      right_contig = pileups.head.rightContig,
      right_start  = pileups.map(_.rightPos).min,
      right_end    = pileups.map(_.rightPos).max,
      right_mean   = Math.round(rightMean).toInt,
      same_strand  = pileups.head.sameStrand,
      mean_count   = meanCount,
      stddev_count = counts.stddev(meanCount),
      median_count = counts.median(),
      pileups      = pileups.length,
      counter      = counter
    )
  }

  val header: IndexedSeq[String] = {
    val names = Metric.names[MergedPileup].filter(_ != "counter")
    val types = EvidenceType.values.map(_.snakeName)
    (names ++ Seq("templates") ++ types).toIndexedSeq
  }

  private val namesAndIndex: Seq[(String, Int)] = {
    Metric.names[MergedPileup]
      .zipWithIndex
      .filter(_._1 != "counter")
      .map { case (name, index) => name -> index }
  }

  def values(mergedPileup: MergedPileup): IndexedSeq[String] = {
    val metricValues = mergedPileup.productIterator.toIndexedSeq
    val staticValues = namesAndIndex.iterator.map { case (_, index) => metricValues(index) }
    val total        = mergedPileup.counter.total
    val counts       = EvidenceType.values.iterator.map { t => mergedPileup.counter.countOf(t) }
    val values       = (staticValues ++ Iterator(total) ++ counts).map(mergedPileup.formatValue).toIndexedSeq
    require(values.length == header.length, s"values: ${values.length} header: ${header.length}")
    values
  }
}