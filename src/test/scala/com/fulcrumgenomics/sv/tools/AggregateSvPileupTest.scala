package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.sv.{BreakpointPileup, UnitSpec}
import com.fulcrumgenomics.util.Metric

import java.nio.file.Path

class AggregateSvPileupTest extends UnitSpec {

  /** Writes the pileups to a temp file and returns the path to the file */
  private def writeTempBreakpointPileups(records: Seq[BreakpointPileup]): Path = {
    val outFile = makeTempFile("pileups", "txt")
    Metric.write[BreakpointPileup](outFile, records)
    outFile
  }

  /** Runs `AggregateSvPileup` end-to-end, reads in the output file and returns its contents as a sequence of records */
  private def runEndToEnd(inputRecords: Seq[BreakpointPileup], maxDist: Int = 100): Seq[AggregatedBreakpointPileup] = {
    val inputFile = writeTempBreakpointPileups(inputRecords)
    val outputFile = makeTempFile("aggregated", "txt")
    new AggregateSvPileup(inputFile, outputFile, maxDist).execute()
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
    split_reads  = 10,
    read_pairs   = 10,
    total        = 20,
  )

  private val pileup_id222_1_100_plus_1_200_plus = BreakpointPileup(
    id           = "222",
    left_contig  = "chr1",
    left_pos     = 100,
    left_strand  = '+',
    right_contig = "chr1",
    right_pos    = 200,
    right_strand = '+',
    split_reads  = 10,
    read_pairs   = 10,
    total        = 20,
  )

  private val pileup_id3_1_100_plus_3_200_minus = BreakpointPileup(
    id           = "3",
    left_contig  = "chr1",
    left_pos     = 100,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 200,
    right_strand = '-',
    split_reads  = 10,
    read_pairs   = 10,
    total        = 20,
  )

  private val pileup_id456_1_200_plus_3_100_plus = BreakpointPileup(
    id           = "456",
    left_contig  = "chr1",
    left_pos     = 200,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 100,
    right_strand = '+',
    split_reads  = 10,
    read_pairs   = 10,
    total        = 20,
  )

  private val pileup_id5_1_300_plus_3_200_plus = BreakpointPileup(
    id           = "5",
    left_contig  = "chr1",
    left_pos     = 300,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 200,
    right_strand = '+',
    split_reads  = 10,
    read_pairs   = 10,
    total        = 20,
  )

  private val pileup_id6_1_300_plus_3_401_plus = BreakpointPileup(
    id           = "6",
    left_contig  = "chr1",
    left_pos     = 300,
    left_strand  = '+',
    right_contig = "chr3",
    right_pos    = 401,
    right_strand = '+',
    split_reads  = 10,
    read_pairs   = 10,
    total        = 20,
  )

  // For example pileups that share contig and strand of both ends, map of pileup to neighbors for use in tests
  private val pileup_to_neighbors_1_plus_plus = Map(
    pileup_id112_1_100_plus_3_200_plus -> Seq(pileup_id456_1_200_plus_3_100_plus),
    pileup_id456_1_200_plus_3_100_plus -> Seq(pileup_id112_1_100_plus_3_200_plus, pileup_id5_1_300_plus_3_200_plus),
    pileup_id5_1_300_plus_3_200_plus -> Seq(pileup_id456_1_200_plus_3_100_plus),
    pileup_id6_1_300_plus_3_401_plus -> Seq(),
  )

  "AggregateSvPileup.toClusters" should "return correct clusters" in {
    val obsClusters = AggregateSvPileup.toClusters(pileup_to_neighbors_1_plus_plus).sortBy(_.size)
    obsClusters.size shouldBe 2
    obsClusters(0) shouldBe Seq(pileup_id6_1_300_plus_3_401_plus)
    obsClusters(1) should contain theSameElementsAs Seq(
      pileup_id112_1_100_plus_3_200_plus, pileup_id456_1_200_plus_3_100_plus, pileup_id5_1_300_plus_3_200_plus
    )
  }

  "AggregateSvPileup.extractClusterFor" should "return correct cluster of size 1" in {
    AggregateSvPileup.extractClusterFor(pileup_id6_1_300_plus_3_401_plus, pileup_to_neighbors_1_plus_plus) should be (
      Seq(pileup_id6_1_300_plus_3_401_plus), pileup_to_neighbors_1_plus_plus.removed(pileup_id6_1_300_plus_3_401_plus)
    )
  }

  it should "return correct cluster of size > 1" in {
    val pileups = Seq(pileup_id112_1_100_plus_3_200_plus, pileup_id456_1_200_plus_3_100_plus, pileup_id5_1_300_plus_3_200_plus)
    pileups.foreach {pileup =>
      AggregateSvPileup.extractClusterFor(pileup, pileup_to_neighbors_1_plus_plus) match {
        case (cluster, _) =>
          cluster should contain theSameElementsAs pileups
      }
    }
  }

  "AggregateSvPileup" should "return an empty seq for empty input" in {
    runEndToEnd(Seq()) shouldEqual Seq()
  }

  it should "aggregate a single pileup" in {
    val pileups = Seq(pileup_id112_1_100_plus_3_200_plus)
    val aggregatedPileups = runEndToEnd(pileups)
    aggregatedPileups should contain theSameElementsAs Seq(
      AggregatedBreakpointPileup(
        id            = "112",
        category      = "Possible inter-contig rearrangement",
        left_contig   = "chr1",
        left_min_pos  = 100,
        left_max_pos  = 100,
        left_strand   = '+',
        right_contig  = "chr3",
        right_min_pos = 200,
        right_max_pos = 200,
        right_strand  = '+',
        split_reads   = 10,
        read_pairs    = 10,
        total         = 20,
      )
    )
  }

  it should "aggregate multiple pileups" in {
    val pileups = Seq(
      pileup_id112_1_100_plus_3_200_plus,
      pileup_id222_1_100_plus_1_200_plus,
      pileup_id3_1_100_plus_3_200_minus,
      pileup_id456_1_200_plus_3_100_plus,
      pileup_id5_1_300_plus_3_200_plus,
      pileup_id6_1_300_plus_3_401_plus,
    )

    // pileup_id112_1_100_plus_3_200_plus, pileup_id456_1_200_plus_3_100_plus, pileup_id5_1_300_plus_3_200_plus
    val expAgg1 = AggregatedBreakpointPileup(
      id            = "112_5_456",
      category      = "Possible inter-contig rearrangement",
      left_contig   = "chr1",
      left_min_pos  = 100,
      left_max_pos  = 300,
      left_strand   = '+',
      right_contig  = "chr3",
      right_min_pos = 100,
      right_max_pos = 200,
      right_strand  = '+',
      split_reads   = 30,
      read_pairs    = 30,
      total         = 60,
    )

    // pileup_id222_1_100_plus_1_200_plus does not combine due to different chromosome
    val expAgg2 = AggregatedBreakpointPileup(
      id            = "222",
      left_contig   = "chr1",
      category      = "Possible deletion",
      left_min_pos  = 100,
      left_max_pos  = 100,
      left_strand   = '+',
      right_contig  = "chr1",
      right_min_pos = 200,
      right_max_pos = 200,
      right_strand  = '+',
      split_reads   = 10,
      read_pairs    = 10,
      total         = 20,
    )

    // pileup_id3_1_100_plus_3_200_minus does not combine due to different strand
    val expAgg3 = AggregatedBreakpointPileup(
      id            = "3",
      category      = "Possible inter-contig rearrangement",
      left_contig   = "chr1",
      left_min_pos  = 100,
      left_max_pos  = 100,
      left_strand   = '+',
      right_contig  = "chr3",
      right_min_pos = 200,
      right_max_pos = 200,
      right_strand  = '-',
      split_reads   = 10,
      read_pairs    = 10,
      total         = 20,
    )

    // pileup_id6_1_300_plus_3_401_plus does not combine due to distance
    val expAgg4 = AggregatedBreakpointPileup(
      id            = "6",
      category      = "Possible inter-contig rearrangement",
      left_contig   = "chr1",
      left_min_pos  = 300,
      left_max_pos  = 300,
      left_strand   = '+',
      right_contig  = "chr3",
      right_min_pos = 401,
      right_max_pos = 401,
      right_strand  = '+',
      split_reads   = 10,
      read_pairs    = 10,
      total         = 20,
    )

    val aggregatedPileups = runEndToEnd(pileups)
    aggregatedPileups should contain theSameElementsAs Seq(expAgg1, expAgg2, expAgg3, expAgg4)

  }

}
