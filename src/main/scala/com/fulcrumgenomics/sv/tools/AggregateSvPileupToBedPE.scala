package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.io.Writer
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sv.tools.AggregateSvPileupToBedPE.BedPE
import com.fulcrumgenomics.sv.tools.AggregateSvPileupToBedPE.BedPE.BedPEWriter
import com.fulcrumgenomics.util.{Io, Metric}
import htsjdk.tribble.annotation.Strand

import java.io.BufferedWriter


@clp(group=ClpGroups.Utilities, description= "Convert the output of AggregateSvPileup to BEDPE.")
class AggregateSvPileupToBedPE(
  @arg(flag='i', doc="Input text file of aggregate pileups generated by AggregateSvPileup") input: FilePath,
  @arg(flag='o', doc="Output text file of the aggregate pileups in BEDPE format.") output: FilePath,
) extends SvTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val reader = Metric.iterator[AggregatedBreakpointPileup](input)
    val writer = BedPEWriter(output)

    reader.map(BedPE.apply).foreach(writer.write)

    writer.close()
  }
}

/** Companion object for [[AggregateSvPileupToBedPE]]. */
object AggregateSvPileupToBedPE {

  /** The IGV-supported BEDPE file extension. */
  val BedPEExtension: String = ".bedpe"

  /** A simple BEDPE record as defined by `bedtools`:
    *
    * - https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
    *
    * Future compatibility could be implemented for supporting 10x flavored BEDPE files:
    *
    * - https://github.com/igvteam/igv/wiki/BedPE-Support
    *
    * Note that the field `score` is allowed to be a string per bedtools!
    */
  case class BedPE(
    chrom1: String,
    start1: Int,
    end1: Int,
    chrom2: String,
    start2: Int,
    end2: Int,
    name: String,
    score: String,
    strand1: Strand,
    strand2: Strand,
  ) extends Metric

  /** Companion object for [[BedPE]]. */
  object BedPE {

    /** Build a [[BedPE]] record from an [[AggregatedBreakpointPileup]]. */
    def apply(pileup: AggregatedBreakpointPileup): BedPE = {
      new BedPE(
        chrom1  = pileup.left_contig,
        start1  = pileup.left_min_pos,
        end1    = pileup.left_max_pos + 1,
        chrom2  = pileup.right_contig,
        start2  = pileup.right_min_pos,
        end2    = pileup.right_max_pos + 1,
        name    = pileup.id,
        score   = pileup.total.toString,
        strand1 = Strand.decode(pileup.left_strand),
        strand2 = Strand.decode(pileup.right_strand),
      )
    }

    /** A writer class for writing [[BedPE]] records. */
    class BedPEWriter(val out: BufferedWriter) extends Writer[BedPE] {

      /** Write a [[BedPE]] record to the underlying writer. */
      override def write(record: BedPE): Unit = {
        out.write(record.chrom1)
        out.write('\t')
        out.write(Integer.toString(record.start1))
        out.write('\t')
        out.write(Integer.toString(record.end1))
        out.write('\t')
        out.write(record.chrom2)
        out.write('\t')
        out.write(Integer.toString(record.start2))
        out.write('\t')
        out.write(Integer.toString(record.end2))
        out.write('\t')
        out.write(record.name)
        out.write('\t')
        out.write(record.score)
        out.write('\t')
        out.write(record.strand1.toString)
        out.write('\t')
        out.write(record.strand2.toString)
        out.newLine()
      }

      /** Closes the underlying writer. */
      override def close(): Unit = out.close()
    }

    /** Companion object to [[BedPEWriter]]. */
    object BedPEWriter {

      /** Constructs a [[BedPEWriter]] that will write to the provided path. */
      def apply(path: PathToIntervals): BedPEWriter = apply(Io.toWriter(path))

      /** Constructs a [[BedPEWriter]] from a [[java.io.Writer]]. */
      def apply(writer: java.io.Writer): BedPEWriter = writer match {
        case bw: BufferedWriter => new BedPEWriter(bw)
        case w                  => new BedPEWriter(new BufferedWriter(w))
      }
    }
  }
}