package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.reflect.ReflectiveBuilder
import com.fulcrumgenomics.commons.util.DelimitedDataParser
import com.fulcrumgenomics.sv.UnitSpec
import com.fulcrumgenomics.sv.tools.AggregateSvPileupToBedPE.BedPE.BedPEWriter
import com.fulcrumgenomics.sv.tools.AggregateSvPileupToBedPE.{BedPE, BedPEExtension}
import com.fulcrumgenomics.util.Metric
import htsjdk.tribble.annotation.Strand


/** Unit tests for [[AggregateSvPileupToBedPE]]. */
class AggregateSvPileupToBedPETest extends UnitSpec {

  /** A test aggregate breakpoint pileup. */
  private val test_aggregate_breakpoint_pileup = AggregatedBreakpointPileup(
    id             = "112",
    category       = "Inter-contig rearrangement",
    left_contig    = "chr1",
    left_min_pos   = 101,
    left_max_pos   = 101,
    left_strand    = '+',
    right_contig   = "chr3",
    right_min_pos  = 201,
    right_max_pos  = 201,
    right_strand   = '-',
    split_reads    = 1,
    read_pairs     = 1,
    total          = 2,
    left_pileups   = PositionList(101),
    right_pileups  = PositionList(201),
  )

  /** A companion test BEDPE record. */
  private val test_bed_pe = BedPE(
    chrom1  = "chr1",
    start1  = 100,
    end1    = 101,
    chrom2  = "chr3",
    start2  = 200,
    end2    = 201,
    name    = "112",
    score   = 2,
    strand1 = Strand.POSITIVE,
    strand2 = Strand.NEGATIVE,
  )

  "AggregateSvPileupToBedPE.BedPE" should "be instantiated from an AggregateBreakpointPileup" in {
    BedPE(test_aggregate_breakpoint_pileup) shouldBe test_bed_pe
  }

  "AggregateSvPileupToBedPE.BedPEWriter" should "write a BedPE record" in {
    val record = new BedPE(
      chrom1  = "chr1",
      start1  = 100,
      end1    = 101,
      chrom2  = "chr3",
      start2  = 200,
      end2    = 201,
      name    = "112",
      score   = 2,
      strand1 = Strand.POSITIVE,
      strand2 = Strand.NEGATIVE,
    )

    val expected = Seq(
      record.chrom1,
      Integer.toString(record.start1),
      Integer.toString(record.end1),
      record.chrom2,
      Integer.toString(record.start2),
      Integer.toString(record.end2),
      record.name,
      record.score.toString,
      record.strand1.toString,
      record.strand2.toString,
    ).toIndexedSeq

    val output = Io.makeTempFile(this.getClass.getSimpleName, BedPEExtension)
    val writer = BedPEWriter(output)
    writer.write(record)
    writer.close()

    val fields: Seq[String] = new ReflectiveBuilder(classOf[BedPE]).argumentLookup.ordered.map(_.name)
    val records = DelimitedDataParser(output, delimiter = '\t', header = fields).toSeq
    records.length shouldBe 1
    val actual = fields.map(field => records.head.get[String](field).value)
    actual should contain theSameElementsInOrderAs expected
  }

  "AggregateSvPileupToBedPE" should "convert an AggregateSvPileup output to a BEDPE file" in {
    val expected = Seq(
      test_bed_pe.chrom1,
      Integer.toString(test_bed_pe.start1),
      Integer.toString(test_bed_pe.end1),
      test_bed_pe.chrom2,
      Integer.toString(test_bed_pe.start2),
      Integer.toString(test_bed_pe.end2),
      test_bed_pe.name,
      test_bed_pe.score.toString,
      test_bed_pe.strand1.toString,
      test_bed_pe.strand2.toString,
    ).toIndexedSeq

    val input  = Io.makeTempFile(this.getClass.getSimpleName, ".txt")
    val output = Io.makeTempFile(this.getClass.getSimpleName, BedPEExtension)
    Metric.write[AggregatedBreakpointPileup](input, test_aggregate_breakpoint_pileup)

    new AggregateSvPileupToBedPE(input = input, output = output).execute()

    val fields: Seq[String] = new ReflectiveBuilder(classOf[BedPE]).argumentLookup.ordered.map(_.name)
    val records = DelimitedDataParser(output, delimiter = '\t', header = fields).toSeq
    records.length shouldBe 1
    val actual = fields.map(field => records.head.get[String](field).value)
    actual should contain theSameElementsInOrderAs expected
  }
}
