package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.PathToBam
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.sv.{BreakpointPileup, UnitSpec}
import com.fulcrumgenomics.testing.SamBuilder
import com.fulcrumgenomics.util.Metric

import java.io.{BufferedWriter, FileWriter}
import java.nio.file.Path

class AggregateSvPileupTest extends UnitSpec {

  /** Writes the pileups to a temp file and returns the path to the file */
  private def writeTempBreakpointPileups(records: Seq[BreakpointPileup]): Path = {
    val outFile = makeTempFile("pileups", "txt")
    Metric.write[BreakpointPileup](outFile, records)
    outFile
  }

  /** Runs `AggregateSvPileup` end-to-end, reads in the output file and returns its contents as a sequence of records */
  private def runEndToEnd(
                           inputRecords: Seq[BreakpointPileup],
                           bam: Option[PathToBam] = None,
                           flank: Int = 500,
                           targets: Option[Path] = None,
                           maxDist: Int = 100): Seq[AggregatedBreakpointPileup] = {
    val inputFile = writeTempBreakpointPileups(inputRecords)
    val outputFile = makeTempFile("aggregated", ".txt")
    new AggregateSvPileup(
      inputFile,
      bam = bam,
      minBreakpointSupport = 0,
      minFrequency = 0,
      flank = flank,
      targetsBed = targets,
      output = outputFile,
      maxDist = maxDist,
    ).execute()
    Metric.read[AggregatedBreakpointPileup](outputFile)
  }

  private val pileup_id112_1_100_plus_3_200_plus = BreakpointPileup(
    id           = "112",
    left_contig  = "chr1",
    left_pos     = 100,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 200,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id222_1_100_plus_1_200_plus = BreakpointPileup(
    id           = "222",
    left_contig  = "chr1",
    left_pos     = 100,
    left_strand  = '+',
    right_contig = "chr1",
    right_pos    = 200,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id3_1_100_plus_3_200_minus = BreakpointPileup(
    id           = "3",
    left_contig  = "chr1",
    left_pos     = 100,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 200,
    right_strand = '-',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id456_1_200_plus_3_100_plus = BreakpointPileup(
    id           = "456",
    left_contig  = "chr1",
    left_pos     = 200,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 100,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id5_1_300_plus_3_200_plus = BreakpointPileup(
    id           = "5",
    left_contig  = "chr1",
    left_pos     = 300,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 200,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id7_1_400_plus_3_300_plus = BreakpointPileup(
    id           = "7",
    left_contig  = "chr1",
    left_pos     = 400,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 300,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id8_1_500_plus_3_400_plus = BreakpointPileup(
    id           = "8",
    left_contig  = "chr1",
    left_pos     = 500,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 400,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id9_1_600_plus_3_500_plus = BreakpointPileup(
    id           = "9",
    left_contig  = "chr1",
    left_pos     = 600,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 500,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id6_1_300_plus_3_401_plus = BreakpointPileup(
    id           = "6",
    left_contig  = "chr1",
    left_pos     = 300,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 401,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val pileup_id10_1_525_plus_3_425_plus = BreakpointPileup(
    id           = "10",
    left_contig  = "chr1",
    left_pos     = 525,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 425,
    right_strand = '+',
    split_reads  = 1,
    read_pairs   = 1,
    total        = 2,
  )

  private val samBuilder = new SamBuilder(readLength = 50, sort = Some(SamOrder.Coordinate),
    sd = Some(SequenceDictionary(
      SequenceMetadata(name = "chr1", index = 0, length = 100000),
      SequenceMetadata(name = "chr2", index = 1, length = 100000),
      SequenceMetadata(name = "chr3", index = 2, length = 100000),
  )))
  samBuilder.addFrag(contig = 0, start = 400)
  samBuilder.addFrag(contig = 0, start = 400)
  samBuilder.addFrag(contig = 0, start = 401)
  samBuilder.addFrag(contig = 0, start = 350)
  samBuilder.addFrag(contig = 0, start = 475)
  samBuilder.addFrag(contig = 0, start = 480)
  samBuilder.addFrag(contig = 0, start = 490)
  samBuilder.addFrag(contig = 0, start = 490)
  samBuilder.addFrag(contig = 0, start = 500)
  samBuilder.addFrag(contig = 0, start = 501)
  samBuilder.addPair(contig = 0, start1 = 320, start2 = 420)
  samBuilder.addPair(contig = 0, start1 = 375, start2 = 401)
  samBuilder.addPair(contig = 0, start1 = 401, start2 = 402)
  samBuilder.addPair(contig = 0, start1 = 480, start2 = 530)
  samBuilder.addPair(contig = 0, start1 = 475, start2 = 476)
  samBuilder.addPair(contig = 0, start1 = 475, start2 = 476)
  samBuilder.addPair(contig = 0, start1 = 500, start2 = 600)
  samBuilder.addPair(contig = 0, start1 = 100, start2 = 700)
  samBuilder.addPair(contig = 0, start1 = 700, start2 = 100)
  samBuilder.addPair(contig = 0, start1 = 100, start2 = 10000)
  samBuilder.addPair(contig = 0, start1 = 10000, start2 = 100)
  samBuilder.addFrag(contig = 2, start = 400)
  samBuilder.addFrag(contig = 2, start = 400)
  samBuilder.addPair(contig = 2, start1 = 375, start2 = 475)
  samBuilder.addPair(contig = 2, start1 = 375, start2 = 475)
  samBuilder.addPair(contig = 2, start1 = 375, start2 = 475)
  samBuilder.addPair(contig = 2, start1 = 375, start2 = 475)
  samBuilder.addPair(contig = 2, start1 = 375, start2 = 475)
  samBuilder.addPair(contig = 2, start1 = 301, start2 = 401)
  private val bam = samBuilder.toTempFile()

  private val bed = makeTempFile("targets", ".bed")
  private val bedWriter = new BufferedWriter(new FileWriter(bed.toFile))
  bedWriter.write("#header\nchr1\t400\t410\nchr3\t401\t402\n")
  bedWriter.flush()
  bedWriter.close()

  // For example pileups that share contig and strand of both ends, map of pileup to neighbors for use in tests
  private val pileup_to_neighbors_1_plus_plus = Map(
    pileup_id112_1_100_plus_3_200_plus -> IndexedSeq(pileup_id456_1_200_plus_3_100_plus),
    pileup_id456_1_200_plus_3_100_plus -> IndexedSeq(pileup_id112_1_100_plus_3_200_plus, pileup_id5_1_300_plus_3_200_plus),
    pileup_id5_1_300_plus_3_200_plus -> IndexedSeq(pileup_id456_1_200_plus_3_100_plus),
    pileup_id6_1_300_plus_3_401_plus -> IndexedSeq(),
  )

  "AggregateSvPileup.toClusters" should "return correct clusters" in {
    val obsClusters = AggregateSvPileup.toClusters(pileup_to_neighbors_1_plus_plus).sortBy(_.size)
    obsClusters.size shouldBe 2
    obsClusters(0) shouldBe IndexedSeq(pileup_id6_1_300_plus_3_401_plus)
    obsClusters(1) should contain theSameElementsAs IndexedSeq(
      pileup_id112_1_100_plus_3_200_plus, pileup_id456_1_200_plus_3_100_plus, pileup_id5_1_300_plus_3_200_plus
    )
  }

  "AggregateSvPileup" should "return an empty seq for empty input" in {
    runEndToEnd(IndexedSeq()) shouldEqual IndexedSeq()
  }

  it should "aggregate a single pileup" in {
    val pileups = IndexedSeq(pileup_id112_1_100_plus_3_200_plus)
    val aggregatedPileups = runEndToEnd(pileups)
    aggregatedPileups should contain theSameElementsAs IndexedSeq(
      AggregatedBreakpointPileup(
        id             = "112",
        category       = "Inter-contig rearrangement",
        left_contig    = "chr1",
        left_min_pos   = 100,
        left_max_pos   = 100,
        left_strand    = '+',
        right_contig   = "chr3",
        right_min_pos  = 200,
        right_max_pos  = 200,
        right_strand   = '+',
        split_reads    = 1,
        read_pairs     = 1,
        total          = 2,
        left_pileups   = PositionList(100),
        right_pileups  = PositionList(200),
      )
    )
  }

  it should "aggregate multiple pileups with no bam or target file" in {
    val pileups = IndexedSeq(
      pileup_id112_1_100_plus_3_200_plus,
      pileup_id222_1_100_plus_1_200_plus,
      pileup_id3_1_100_plus_3_200_minus,
      pileup_id456_1_200_plus_3_100_plus,
      pileup_id5_1_300_plus_3_200_plus,
      pileup_id7_1_400_plus_3_300_plus,
      pileup_id8_1_500_plus_3_400_plus,
      pileup_id9_1_600_plus_3_500_plus,
      pileup_id6_1_300_plus_3_401_plus,
    )

    val expAgg1 = AggregatedBreakpointPileup(
      id             = "112_456_5_7_8_9",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 100,
      left_max_pos   = 600,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 100,
      right_max_pos  = 500,
      right_strand   = '+',
      split_reads    = 6,
      read_pairs     = 6,
      total          = 12,
      left_pileups   = PositionList(100, 200, 300, 400, 500, 600),
      right_pileups  = PositionList(100, 200, 200, 300, 400, 500),
    )

    // pileup_id222_1_100_plus_1_200_plus does not combine due to different chromosome
    val expAgg2 = AggregatedBreakpointPileup(
      id             = "222",
      left_contig    = "chr1",
      category       = "Possible deletion",
      left_min_pos   = 100,
      left_max_pos   = 100,
      left_strand    = '+',
      right_contig   = "chr1",
      right_min_pos  = 200,
      right_max_pos  = 200,
      right_strand   = '+',
      split_reads    = 1,
      read_pairs     = 1,
      total          = 2,
      left_pileups   = PositionList(100),
      right_pileups  = PositionList(200),
    )

    // pileup_id3_1_100_plus_3_200_minus does not combine due to different strand
    val expAgg3 = AggregatedBreakpointPileup(
      id             = "3",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 100,
      left_max_pos   = 100,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 200,
      right_max_pos  = 200,
      right_strand   = '-',
      split_reads    = 1,
      read_pairs     = 1,
      total          = 2,
      left_pileups   = PositionList(100),
      right_pileups  = PositionList(200),
    )

    // pileup_id6_1_300_plus_3_401_plus does not combine due to distance
    val expAgg4 = AggregatedBreakpointPileup(
      id             = "6",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 300,
      left_max_pos   = 300,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 401,
      right_max_pos  = 401,
      right_strand   = '+',
      split_reads    = 1,
      read_pairs     = 1,
      total          = 2,
      left_pileups   = PositionList(300),
      right_pileups  = PositionList(401),
    )

    val aggregatedPileups = runEndToEnd(pileups)
    aggregatedPileups should contain theSameElementsAs IndexedSeq(expAgg1, expAgg2, expAgg3, expAgg4)
  }

  it should "aggregate multiple pileups with bam and target files" in {
    val pileups = IndexedSeq(
      pileup_id6_1_300_plus_3_401_plus,
      pileup_id7_1_400_plus_3_300_plus,
      pileup_id8_1_500_plus_3_400_plus,
      pileup_id9_1_600_plus_3_500_plus,
      pileup_id10_1_525_plus_3_425_plus,
    )

    val expAgg1 = AggregatedBreakpointPileup(
      id             = "10_7_8_9",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 400,
      left_max_pos   = 600,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 300,
      right_max_pos  = 500,
      right_strand   = '+',
      split_reads    = 4,
      read_pairs     = 4,
      total          = 8,
      left_pileups   = PositionList(400, 500, 525, 600),
      right_pileups  = PositionList(300, 400, 425, 500),
      left_frequency = Some(8d/16),
      right_frequency = Some(8d/8),
      left_overlaps_target = true,
      right_overlaps_target = true,
    )

    val expAgg2 = AggregatedBreakpointPileup(
      id             = "6",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 300,
      left_max_pos   = 300,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 401,
      right_max_pos  = 401,
      right_strand   = '+',
      split_reads    = 1,
      read_pairs     = 1,
      total          = 2,
      left_pileups   = PositionList(300),
      right_pileups  = PositionList(401),
      left_frequency = Some(2d/2),
      right_frequency = Some(2d/8),
      left_overlaps_target = false,
      right_overlaps_target = false,  // There is a BED feature 401-402; BED is 0-based half open
    )

    val aggregatedPileups = runEndToEnd(pileups, bam = Some(bam), targets = Some(bed))
    aggregatedPileups should contain theSameElementsAs IndexedSeq(expAgg1, expAgg2)

  }

  it should "aggregate multiple pileups with a bam file" in {
    val pileups = IndexedSeq(
      pileup_id7_1_400_plus_3_300_plus,
      pileup_id8_1_500_plus_3_400_plus,
      pileup_id9_1_600_plus_3_500_plus,
      pileup_id10_1_525_plus_3_425_plus,
    )


    val expAgg = AggregatedBreakpointPileup(
      id             = "10_7_8_9",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 400,
      left_max_pos   = 600,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 300,
      right_max_pos  = 500,
      right_strand   = '+',
      split_reads    = 4,
      read_pairs     = 4,
      total          = 8,
      left_pileups   = PositionList(400, 500, 525, 600),
      right_pileups  = PositionList(300, 400, 425, 500),
      left_frequency = Some(8d/16),
      right_frequency = Some(8d/8),
    )

    val aggregatedPileups = runEndToEnd(pileups, bam = Some(bam))
    aggregatedPileups should contain theSameElementsAs IndexedSeq(expAgg)

  }

  it should "aggregate multiple pileups with a target file" in {
    val pileups = IndexedSeq(
      pileup_id7_1_400_plus_3_300_plus,
      pileup_id8_1_500_plus_3_400_plus,
      pileup_id9_1_600_plus_3_500_plus,
      pileup_id10_1_525_plus_3_425_plus,
    )

    val expAgg = AggregatedBreakpointPileup(
      id             = "10_7_8_9",
      category       = "Inter-contig rearrangement",
      left_contig    = "chr1",
      left_min_pos   = 400,
      left_max_pos   = 600,
      left_strand    = '+',
      right_contig   = "chr3",
      right_min_pos  = 300,
      right_max_pos  = 500,
      right_strand   = '+',
      split_reads    = 4,
      read_pairs     = 4,
      total          = 8,
      left_pileups   = PositionList(400, 500, 525, 600),
      right_pileups  = PositionList(300, 400, 425, 500),
      left_frequency = None,
      right_frequency = None,
      left_overlaps_target = true,
      right_overlaps_target = true,
    )

    val aggregatedPileups = runEndToEnd(pileups, targets = Some(bed))
    aggregatedPileups should contain theSameElementsAs IndexedSeq(expAgg)

  }

}
