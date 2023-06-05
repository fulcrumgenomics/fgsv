package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam}
import com.fulcrumgenomics.bam.api.{QueryType, SamRecord, SamSource}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sv.tools.AggregateSvPileup.PileupGroup
import com.fulcrumgenomics.sv.tools.BreakpointCategory.BreakpointCategory
import com.fulcrumgenomics.sv.{BreakpointPileup, Intervals}
import com.fulcrumgenomics.util.{Io, Metric}
import htsjdk.samtools.util.{Interval, OverlapDetector}
import htsjdk.tribble.bed.BEDFeature
import com.fulcrumgenomics.FgBioDef.javaIterableToIterator
import scala.collection.mutable


@clp(group=ClpGroups.All, description=
  """
    |Merges nearby pileups of reads supporting putative breakpoints.
    |
    |Takes as input the file of pileups produced by `SvPileup`. That file contains a list of breakpoints, each
    |consisting of a chromosome, position and strand for each side of the breakpoint, as well as quantified read support
    |for the breakpoint.
    |
    |This tool merges sets of breakpoints that have their left sides on the same chromosome, their right sides on the
    |same chromosome, the same left and right strands, and their left and right positions both within a length
    |threshold. The merging behavior is transitive. For example, two breakpoints that are farther apart than the length
    |threshold can be merged if there is a third breakpoint that is close to both.
    |
    |`SvPileup` distinguishes between two types of evidence for breakpoint events: split-read evidence, where a
    |breakpoint occurs within a mapped read, and read-pair evidence, where a breakpoint occurs in the unsequenced
    |insert between two paired reads. Currently this tool treats both forms of evidence equally, despite the
    |inaccuracy of positions reported by `SvPileup` for read-pair evidence.
    |
    |If a BAM file is provided, each aggregated pileup is annotated with the allele frequency at its left and right
    |breakends. Allele frequency is defined as the total number of templates supporting constituent breakpoints divided
    |by the count of the union of templates that cross any constituent breakends. In particular, paired templates that
    |straddle a breakend are considered to cross the breakend.
    |
    |If a BED file of target regions is provided, each aggregated pileup is annotated with whether its left and
    |right sides overlap a target region.  Additionally, the names of the overlapping target regions will be copied
    |annotated, overwriting any values from the `SvPileup` input.  If no target regions are provided, then the names
    |of the overlapping target regions are copied from the `SvPiluep` input (if present).
    |
    |The output file is a tab-delimited table with one record per aggregated cluster of pileups. Aggregated
    |pileups are reported with the minimum and maximum (inclusive) coordinates of all pileups in the cluster, a
    |possible putative structural variant event type supported by the pileups, and the sum of read support from all
    |pileups in the cluster.
    |""")
class AggregateSvPileup
(@arg(flag='i', doc="Input text file of pileups generated by SvPileup") input: FilePath,
 @arg(flag='b', doc="Bam file for allele frequency analysis. Must be the same file that was used as input to SvPileup.")
 bam: Option[PathToBam] = None,
 @arg(flag='f', doc="If BAM file is provided: distance upstream and downstream of aggregated breakpoint regions to " +
   "search for mapped templates that overlap breakends. These are the templates that will be partitioned into those " +
   "supporting the breakpoint vs. reading through it for the allele frequency calculation. Recommended to use at " +
   "least the max read pair inner distance used by SvPileup.") flank: Int = 1000,
 @arg(doc="If BAM file is provided: minimum total number of templates supporting an aggregated breakpoint to report " +
   "allele frequency. Supports speed improvement by avoiding querying and iterating over huge read pileups that " +
   "contain insufficient support for a breakpoint to be considered interesting.") minBreakpointSupport: Int = 10,
 @arg(doc="If BAM file is provided: minimum allele frequency to report. Supports speed improvement by " +
   "avoiding iterating over huge read pileups that contain insufficient support for a breakpoint to be considered " +
   "interesting.") minFrequency: Double = 0.001,
 @arg(flag='t', doc="Optional BED file of target regions") targetsBed: Option[FilePath] = None,
 @arg(flag='o', doc="Output file") output: FilePath,
 @arg(flag='d', doc="Distance threshold below which to cluster breakpoints") maxDist: Int = 10,
) extends SvTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {

    // Read pileups from input file
    val pileups: PileupGroup = Metric.read[BreakpointPileup](input).toIndexedSeq
    logger.info(f"Read ${pileups.length} pileups from file $input")

    // Open BAM reader
    val samSource: Option[SamSource] = bam.map(SamSource(_))

    // Read targets from BED file and build OverlapDetector
    val targets: Option[OverlapDetector[BEDFeature]] = targetsBed.map(Intervals.overlapDetectorFrom)

    // Open output writer
    val writer = Metric.writer[AggregatedBreakpointPileup](output)

    // Group pileups by left and right contig and strand
    pileups.groupBy(_.contigOrientation).foreach { case (_, group) =>

      // Within each contig+strand group, identify all immediate "neighbors" of each pileup (other pileups within the
      // distance threshold)
      val pileupToNeighbors: Map[BreakpointPileup, PileupGroup] = AggregateSvPileup.pileupToNeighbors(group, maxDist)

      // Within each contig+strand group, recursively build (transitive) clusters of neighboring pileups
      val pileupClusters: Seq[PileupGroup] = AggregateSvPileup.toClusters(pileupToNeighbors)
      logger.info(f"Processing ${pileupClusters.size} pileup clusters with endpoints " +
        f"${group.head.left_contig}:${group.head.left_strand}, ${group.head.right_contig}:${group.head.right_strand}")

      // Aggregate each cluster, analyze allele frequency and target overlap, and write output
      pileupClusters
        .map(AggregatedBreakpointPileup.apply)
        .map(_.calculateFrequency(samSource, flank, minBreakpointSupport, minFrequency))
        .map(_.addTargetOverlap(targets))
        .foreach(writer.write)
    }

    writer.close()
  }
}


object AggregateSvPileup {

  /** Type alias for a collection of pileups */
  type PileupGroup = IndexedSeq[BreakpointPileup]

  /**
   * Maps each pileup to the set of other pileups within the distance threshold.
   *
   * Note: it is recommended to call this method on a group of pileups that share a common `BreakpointContigOrientation`
   * to reduce unnecessary run time.
   *
   * @param pileups Set of pileups to evaluate for pairwise neighbor relationships
   * @param maxDist Maximum distance (inclusive) for the left and right coordinates to be considered a neighbor
   * @return Map of pileup to the set of neighbors
   */
  private def pileupToNeighbors(pileups: PileupGroup, maxDist: Int): Map[BreakpointPileup, PileupGroup] = {
    val map: mutable.Map[BreakpointPileup, PileupGroup] = mutable.Map.from(pileups.map((_, IndexedSeq())))
    for (i <- pileups.indices; j <- pileups.indices.drop(i)) {
      val pileup1 = pileups(i)
      val pileup2 = pileups(j)
      if (isCluster(pileup1, pileup2, maxDist)) {
        map.put(pileup1, pileup2 +: map(pileup1))
        map.put(pileup2, pileup1 +: map(pileup2))
      }
    }
    map.toMap
  }

  /** Returns true if the two breakpoints can be clustered, i.e., they are not the same breakpoint (have different
   * IDs), they have the same left and right chromosomes and strand, and the distance between their left and right
   * coordinates is below the threshold (inclusive) */
  private def isCluster(bp1: BreakpointPileup, bp2: BreakpointPileup, maxDist: Int): Boolean = {
    bp1.id != bp2.id &&
      bp1.contigOrientation == bp2.contigOrientation &&
      Math.abs(bp1.left_pos - bp2.left_pos) <= maxDist &&
      Math.abs(bp1.right_pos - bp2.right_pos) <= maxDist
  }

  /**
   * Converts a set of pairwise neighbor relationships to a set of clusters, where each cluster consists of breakpoints
   * that can be joined into a connected path via pairwise neighbor relationships.
   * @param pileupToNeighbors Map of pileup to its set of neighbors
   * @return The set of pileup clusters
   */
  def toClusters(pileupToNeighbors: Map[BreakpointPileup, PileupGroup]): Seq[PileupGroup] = {
    // The following finds connnected components in a breadth first manner
    val visited    = scala.collection.mutable.HashSet[BreakpointPileup]()
    val components = IndexedSeq.newBuilder[PileupGroup]
    pileupToNeighbors.keys.foreach { rootPileup =>
      if (!visited.contains(rootPileup)) {
        val remaining = scala.collection.mutable.Queue[BreakpointPileup]()
        val component = Set.newBuilder[BreakpointPileup]
        remaining.enqueue(rootPileup)
        while (remaining.nonEmpty) {
          // get the next pileup to examine
          val curPileup = remaining.dequeue()
          // add it to this component
          component += curPileup
          // find all unvisited neighbors
          val unvisited = pileupToNeighbors(curPileup).filterNot(visited.contains)
          // add the unvisited neighbors to the queue
          remaining ++= unvisited
          // mark the unvisited neighbors as being visited so we don't enqueue them later
          unvisited.foreach(visited.add)
          // mark this pileup as visited
          visited.add(curPileup)
        }
        components += component.result().toIndexedSeq
      }
    }

    components.result()
  }
}


/** Contig and strand for the left and right sides of a breakpoint pileup */
case class BreakpointContigOrientation(left_contig: String, left_strand: Char, right_contig: String, right_strand: Char)


object BreakpointContigOrientation {
  def apply(pileup: BreakpointPileup): BreakpointContigOrientation = {
    BreakpointContigOrientation(
      left_contig  = pileup.left_contig,
      left_strand  = pileup.left_strand,
      right_contig = pileup.right_contig,
      right_strand = pileup.right_strand,
    )
  }
}


/** Type of structural variant supported by a pileup */
object BreakpointCategory extends Enumeration {
  type BreakpointCategory = String
  val PossibleDeletion = "Possible deletion"
  val IntraContig = "Intra-contig rearrangement"
  val InterContig = "Inter-contig rearrangement"

  def apply(breakpointPileup: BreakpointPileup): BreakpointCategory = {
    if (breakpointPileup.left_contig != breakpointPileup.right_contig) {
      InterContig
    } else {
      if (breakpointPileup.left_strand != breakpointPileup.right_strand) {
        IntraContig
      } else PossibleDeletion
    }
  }
}


/**
 * Aggregated cluster of breakpoint pileups
 * @param id                    Combined ID retaining the IDs of all constituent breakpoints
 * @param category              Breakpoint category
 * @param left_contig           Contig name for left side of breakpoint
 * @param left_min_pos          Minimum coordinate of left breakends
 * @param left_max_pos          Maximum coordinate of left breakends
 * @param left_strand           Strand at left breakends
 * @param right_contig          Contig name for right side of breakpoint
 * @param right_min_pos         Minimum coordinate of right breakends
 * @param right_max_pos         Maximum coordinate of right breakends
 * @param right_strand          Strand at right breakends
 * @param split_reads           Total number of split reads supporting the breakpoints in the cluster
 * @param read_pairs            Total number of read pairs supporting the breakpoints in the cluster
 * @param total                 Total number of templates supporting the breakpoints in the cluster
 * @param left_pileups          List of constituent left breakends
 * @param right_pileups         List of constituent right breakends
 * @param left_frequency        Proportion of templates mapping at one of the left breakends that support a breakpoint.
 *                              If a template maps across multiple breakends, it is only counted once. If a template
 *                              maps entirely between two breakends but does not overlap one, it is not counted. If two
 *                              paired reads straddle a breakend, the template is counted as reading across the
 *                              breakend. If two paired reads both overlap the same breakend, the template is counted
 *                              once.
 * @param right_frequency       Proportion of reads mapping at one of the right breakends that support a breakpoint
 * @param left_overlaps_target  True if the left aggregated region overlaps a target region
 * @param right_overlaps_target True if the right aggregated region overlaps a target region
 * @param left_targets  the comma-delimited list of target names overlapping the left breakpoint
 * @param right_targets the comma-delimited list of target names overlapping the right breakpoint
 */
case class AggregatedBreakpointPileup(id: String,
                                      category: BreakpointCategory,
                                      left_contig: String,
                                      left_min_pos: Int,
                                      left_max_pos: Int,
                                      left_strand: Char,
                                      right_contig: String,
                                      right_min_pos: Int,
                                      right_max_pos: Int,
                                      right_strand: Char,
                                      split_reads: Int,
                                      read_pairs: Int,
                                      total: Int,
                                      left_pileups: PositionList,
                                      right_pileups: PositionList,
                                      left_frequency: Option[Double] = None,
                                      right_frequency: Option[Double] = None,
                                      left_overlaps_target: Boolean = false,
                                      right_overlaps_target: Boolean = false,
                                      left_targets: Option[String] = None,
                                      right_targets: Option[String] = None
                                     ) extends Metric with LazyLogging {

  /** Returns the IDs of constituent breakpoints */
  def pileupIds(): Seq[String] = id.split("_").toSeq.sorted

  /** Returns a new aggregated pileup with the given pileup added */
  def add(pileup: BreakpointPileup): AggregatedBreakpointPileup = {
    assert(pileup.left_contig == left_contig)
    assert(pileup.right_contig == right_contig)
    assert(pileup.left_strand == left_strand)
    assert(pileup.right_strand == right_strand)

    val left_targets = (this.left_targets, pileup.left_targets) match {
      case (_, None) => this.left_targets
      case (None, _) => pileup.left_targets
      case (Some(this_left), Some(pileup_left)) => Some(f"${this_left},${pileup_left}")
    }
    val right_targets = (this.right_targets, pileup.right_targets) match {
      case (_, None) => this.right_targets
      case (None, _) => pileup.right_targets
      case (Some(this_right), Some(pileup_right)) => Some(f"${this_right},${pileup_right}")
    }

    AggregatedBreakpointPileup(
      id              = pileupIds().appended(pileup.id).sorted.mkString("_"),
      category        = category,
      left_contig     = left_contig,
      left_min_pos    = Math.min(left_min_pos, pileup.left_pos),
      left_max_pos    = Math.max(left_max_pos, pileup.left_pos),
      left_strand     = left_strand,
      right_contig    = right_contig,
      right_min_pos   = Math.min(right_min_pos, pileup.right_pos),
      right_max_pos   = Math.max(right_max_pos, pileup.right_pos),
      right_strand    = right_strand,
      split_reads     = split_reads + pileup.split_reads,
      read_pairs      = read_pairs + pileup.read_pairs,
      total           = total + pileup.total,
      left_pileups    = left_pileups.appended(pileup.left_pos),
      right_pileups   = right_pileups.appended(pileup.right_pos),
      left_targets    = left_targets,
      right_targets   = right_targets,
    )
  }

  /**
   * Returns a new aggregated pileup with left and right frequency populated
   * @param source Optional BAM reader; frequencies are set to None if not provided
   * @param flank Distance upstream and downstream of this aggregated pileup region to query for templates that read
   *              across breakends
   * @param minTotal Minimum total number of templates supporting the aggregated breakpoint in order to calculate
   *                 frequency. Used to speed up execution when the number of supporting templates is too low to be
   *                 interesting.
   * @param minFrequency Minimum proportion of templates overlapping a breakend that support the breakpoint. Used to
   *                     stop iteration over huge sets of overlapping templates when the number of templates supporting
   *                     the breakpoint is insufficient to be interesting. When the number of overlapping templates
   *                     causes (total_supporting_breakpoint / total_overlappers) to drop below this frequency, None is
   *                     returned.
   */
  def calculateFrequency(source: Option[SamSource],
                         flank: Int,
                         minTotal: Int,
                         minFrequency: Double):
  AggregatedBreakpointPileup = {

    def frequency(contig: String, pileups: PositionList, minPos: Int, maxPos: Int): Option[Double] = source match {
      case None => None
      case Some(src) =>
        val numOverlap = numOverlappingTemplates(
          source       = src,
          contig       = contig,
          breakends    = pileups,
          queryMin     = minPos - flank,
          queryMax     = maxPos + flank,
          minTotal     = minTotal,
          minFrequency = minFrequency,
        )
        numOverlap match {
          case None => None
          case Some(n) => Some(total.toDouble / n)
        }
    }

    this.copy(
      left_frequency = frequency(left_contig, left_pileups, left_min_pos, left_max_pos),
      right_frequency = frequency(right_contig, right_pileups, right_min_pos, right_max_pos),
    )
  }

  /**
   * Returns the number of templates that overlap any breakend in the list. Returns None if there are insufficient
   * templates supporting the breakpoint, or if the proportion of those templates that support the breakpoint is less
   * than the minimum frequency.
   * @param source BAM reader
   * @param contig Contig name
   * @param breakends Breakend positions with coordinates defined as in `SvPileup` (1-based)
   * @param queryMin Minimum coordinate to query for overlapping templates (1-based)
   * @param queryMax Maximum coordinate to query for overlapping templates (1-based)
   * @param minTotal Minimum total number of templates supporting the aggregated breakpoint. Used to speed up
   *                 execution when the number of supporting templates is too low to be interesting.
   * @param minFrequency Minimum proportion of templates overlapping a breakend that support the breakpoint. Used to
   *                     stop iteration over huge sets of overlapping templates when the number of templates supporting
   *                     the breakpoint is insufficient to be interesting. When the number of overlapping templates
   *                     causes (total_supporting_breakpoint / total_overlappers) to drop below this frequency, None is
   *                     returned.
   * @return The size of the union of templates that overlap any breakend
   */
  private def numOverlappingTemplates(source: SamSource,
                                      contig: String,
                                      breakends: PositionList,
                                      queryMin: Int,
                                      queryMax: Int,
                                      minTotal: Int,
                                      minFrequency: Double): Option[Int] = {

    if (total < minTotal) return None

    // Read names of overlapping templates
    val overlappers: mutable.Set[String] = new mutable.HashSet[String]()

    // Max number of overlappers; when this number is exceeded (causing the breakpoint frequency to drop below the
    // minimum), None is returned
    val maxOverlappers = total.toDouble / minFrequency

    // Whether to check the full template for overlap including both reads and insert
    def checkAsPair(rec: SamRecord): Boolean = {
      rec.paired && rec.mapped && rec.mateMapped &&
        rec.refName == rec.mateRefName &&
        ((rec.start <= rec.mateStart && rec.positiveStrand && rec.mateNegativeStrand)
          || (rec.start >= rec.mateStart && rec.negativeStrand && rec.matePositiveStrand))
    }

    // Iterate over all reads overlapping the region
    val samIter = source.query(contig, queryMin, queryMax, QueryType.Overlapping)
    while (samIter.hasNext) {
      if (overlappers.size > maxOverlappers) {
        samIter.close()
        return None
      }
      val rec = samIter.next()
      if (!overlappers.contains(rec.name)) {
        if (checkAsPair(rec)) {
          if (breakends.positions.exists(
            pos => pos >= Math.min(rec.start, rec.mateStart) && pos <= Math.max(rec.end, rec.mateEnd.get))) {
            overlappers.add(rec.name)
          }
        } else if (breakends.positions.exists(pos => pos >= rec.start && pos <= rec.end)) {
          overlappers.add(rec.name)
        }
      }
    }

    Some(overlappers.size)

  }

  /**
   * Returns a new aggregated pileup with target overlap fields populated
   * @param targets Optional OverlapDetector for target intervals. If None, target overlap fields are set to false.
   */
  def addTargetOverlap(targets: Option[OverlapDetector[BEDFeature]]): AggregatedBreakpointPileup = {

    def annotations(contig: String, minPos: Int , maxPos: Int): (Boolean, Option[String]) = {
      targets match {
        case None => (false, None)
        case Some(detector) => {
          val span = new Interval(contig, minPos, maxPos)
          if (detector.overlapsAny(span)) {
            val overlaps = detector.getOverlaps(span)
            (true, Some(overlaps.map(_.getName).toSeq.sorted.distinct.mkString(",")))
          } else {
            (false, None)
          }
        }
      }
    }

    val (leftOverlaps, leftTargets) = annotations(left_contig, left_min_pos, left_max_pos)
    val (rightOverlaps, rightTargets) = annotations(right_contig, right_min_pos, right_max_pos)
    this.copy(
      left_overlaps_target  = leftOverlaps,
      right_overlaps_target = rightOverlaps,
      left_targets          = leftTargets,
      right_targets         = rightTargets
    )
  }

  private def overlapsTarget(contig: String, start: Int, end: Int, targets: OverlapDetector[BEDFeature]): Boolean = {
    targets.overlapsAny(new Interval(contig, start, end))
  }

}


object AggregatedBreakpointPileup {

  def apply(pileups: PileupGroup): AggregatedBreakpointPileup = {
    pileups match {
      case head +: tail =>
        val headAgg = AggregatedBreakpointPileup(
          id             = head.id,
          category       = BreakpointCategory(head),
          left_contig    = head.left_contig,
          left_min_pos   = head.left_pos,
          left_max_pos   = head.left_pos,
          left_strand    = head.left_strand,
          right_contig   = head.right_contig,
          right_min_pos  = head.right_pos,
          right_max_pos  = head.right_pos,
          right_strand   = head.right_strand,
          split_reads    = head.split_reads,
          read_pairs     = head.read_pairs,
          total          = head.total,
          left_pileups   = PositionList(head.left_pos),
          right_pileups  = PositionList(head.right_pos),
          left_targets   = head.left_targets,
          right_targets  = head.right_targets,
        )
        tail.foldLeft[AggregatedBreakpointPileup](headAgg)((agg, nextPileup) => agg.add(nextPileup))
      case _ => throw new IllegalArgumentException("Pileup group must be non-empty")
    }
  }

}


/** Thin wrapper around a Seq[Int] primarily used to convert back and forth to a string representation */
class PositionList(val positions: Seq[Int]) {
  override def equals(o: Any): Boolean = o match {
    case p: PositionList => p.positions == this.positions
    case _ => false
  }

  override def hashCode(): Int = positions.hashCode()

  override def toString(): String = positions.mkString(",")

  def appended(pos: Int): PositionList = new PositionList(positions.appended(pos).sorted)
}

object PositionList {
  def apply(positions: Int*): PositionList = new PositionList(positions)
  def apply(s: String): PositionList = new PositionList(s.split(",").map(_.toInt))
}
