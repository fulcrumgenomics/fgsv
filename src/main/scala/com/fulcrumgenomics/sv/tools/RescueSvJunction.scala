package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam}
import com.fulcrumgenomics.bam.api.{QueryType, SamSource}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sv.{Breakpoint, JunctionRescuer}
import com.fulcrumgenomics.util.{Io, Metric}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools.reference.{FastaSequenceIndex, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.Interval

import java.nio.file.{Files, Path}
import java.nio.file.StandardCopyOption.{ATOMIC_MOVE, REPLACE_EXISTING}
import scala.collection.mutable
import scala.util.Using

/** Experimental candidate-guided realignment, deliberately separate from pileup counts. */
@clp(group=ClpGroups.BreakpointAndSv, description="""
  |Experimental soft-clip rescue for one supplied, oriented junction. Compares a read against
  |the synthetic junction and both local reference alternatives, requiring aligned anchors on
  |both sides. Reports qualifying realigned support, including existing discovery supporters;
  |consumers must union template identities with original evidence, not add counts. Does not
  |modify the BAM, emit SA tags, discover junctions, or update pileup counts or frequencies.
  |
  |Eligible reads must be mapped primary alignments, not secondary, supplementary, duplicate,
  |or vendor-QC-failed, and meet --min-mapq. A terminal soft clip must face the supplied endpoint
  |within --slop. Both the observed clip and the immediately adjacent contiguous aligned source
  |segment must meet --min-anchor; the final junction alignment must also meet --min-anchor on
  |each side. Hard clipping does not count as observed support. Equivalent adjacent CIGAR
  |operations are coalesced for eligibility, and padding does not interrupt a source anchor.
  |SAM SEQ '=' is expanded from the original aligned reference positions; '=' in soft clips
  |or insertions has no reference base and must be replaced with literal bases before rescue.
  |
  |Output is a headered metrics table, including when no templates qualify. Each row identifies
  |the supplied endpoints and the actual synthetic-junction boundary after genomic-edge clamping.
  |Coordinates are 1-based inclusive. Query coordinates use the orientation in query_reversed;
  |target coordinates refer to the oriented synthetic junction. Endpoint strand conventions
  |match SvPileup's left/from and right/into endpoints. For a negative strand use --left-strand=-
  |or --right-strand=-. Deduplication applies only within this supplied junction.
  |For example, with flank 200, --left-position 600 --left-strand=+ --right-position 401
  |--right-strand=+ joins reference bases 401..600 to 401..600 at target boundary 200.
  |With --left-strand=- --right-strand=- at the same endpoints, it joins reverse-complemented
  |600..799 to reverse-complemented 202..401, also at boundary 200. Contigs are supplied separately.
  |
  |The representative per template is selected by descending score advantage, junction score,
  |and minimum anchor, then ascending read number, genomic/alignment coordinates, CIGARs and
  |query orientation (false first), with remaining diagnostics breaking ties. This does not
  |override rejection of ambiguous optimal alignment placements.
  |
  |Smith-Waterman scores are fixed: canonical match +2, mismatch/noncanonical base -3,
  |gap open -6 and gap extension -1. Anchors count canonical aligned columns, including
  |mismatches; gaps and ambiguous bases do not contribute. Scores compare only the supplied
  |local alternatives, not all genomic mappings, and are not quality-aware genotype likelihoods
  |or supplementary mapping qualities. The flank must cover each otherwise eligible read's
  |observed length plus slop; increase --flank for longer reads. Genomic edges are clamped,
  |not compensated by truncating reads.
  |""")
class RescueSvJunction
(@arg(flag='i', doc="Indexed coordinate-sorted original BAM") input: PathToBam,
 @arg(flag='r', doc="Indexed reference FASTA used for the BAM") reference: FilePath,
 @arg(flag='o', doc="Output qualifying-template metrics; existing discovery supporters are included") output: FilePath,
 @arg(doc="Left/from reference contig") leftContig: String,
 @arg(doc="Left/from aligned endpoint, 1-based inclusive") leftPosition: Int,
 @arg(doc="Left/from strand, + or -") leftStrand: String,
 @arg(doc="Right/into reference contig") rightContig: String,
 @arg(doc="Right/into aligned endpoint, 1-based inclusive") rightPosition: Int,
 @arg(doc="Right/into strand, + or -") rightStrand: String,
 @arg(doc="Minimum aligned bases on each side of the junction") minAnchor: Int = 10,
 @arg(doc="Maximum distance from a soft-clip boundary to its candidate endpoint") slop: Int = 5,
 @arg(doc="Minimum primary alignment mapping quality") minMapq: Int = 30,
 @arg(doc="Minimum junction score advantage over both local reference alternatives") minScoreDelta: Int = 10,
 @arg(doc="Minimum junction-alignment identity") minIdentity: Double = 0.9,
 @arg(doc="Minimum fraction of observed query bases included in the junction alignment") minQueryCoverage: Double = 0.9,
 @arg(doc="Reference flank and BAM query padding; must cover eligible read length plus slop") flank: Int = 250
) extends SvTool {
  Io.assertReadable(input)
  Io.assertReadable(reference)
  assertOutputSafe()
  Io.assertCanWriteFile(output)
  require(Set("+", "-").contains(leftStrand), "--left-strand must be + or -")
  require(Set("+", "-").contains(rightStrand), "--right-strand must be + or -")
  require(minAnchor > 0, "--min-anchor must be positive")
  require(slop >= 0, "--slop must be non-negative")
  require(flank >= minAnchor && flank > slop, "Increase --flank: it must be at least --min-anchor and greater than --slop")
  require(minMapq >= 0 && minMapq <= 255, "--min-mapq must be in [0,255]")
  require(minScoreDelta >= 0, "--min-score-delta must be non-negative")
  require(minIdentity > 0 && minIdentity <= 1, "--min-identity must be in (0,1]")
  require(minQueryCoverage > 0 && minQueryCoverage <= 1, "--min-query-coverage must be in (0,1]")

  /** Protect inputs and their sidecars, including symlink and hardlink output aliases. */
  private def assertOutputSafe(): Unit = {
    require(!Files.isSymbolicLink(output) || Files.exists(output),
      s"--output must not be a dangling symbolic link ($output)")

    def appended(path: Path, suffix: String): Path = path.resolveSibling(path.getFileName.toString + suffix)

    def bamSidecars(path: Path): IndexedSeq[Path] = {
      val basename = path.resolveSibling(path.getFileName.toString.stripSuffix(".bam"))
      IndexedSeq(".bai", ".csi").flatMap(suffix => IndexedSeq(appended(path, suffix), appended(basename, suffix)))
    }

    def referenceSidecars(path: Path): IndexedSeq[Path] = IndexedSeq(
      ReferenceSequenceFileFactory.getFastaIndexFileName(path),
      appended(path, ".gzi"),
      ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(path)
    )

    val bamPaths = IndexedSeq(input, input.toRealPath()).distinct
    val referencePaths = IndexedSeq(reference, reference.toRealPath()).distinct
    val protectedPaths = IndexedSeq(input -> "--input", reference -> "--reference") ++
      bamPaths.flatMap(bamSidecars).map(_ -> "--input index") ++
      referencePaths.flatMap(referenceSidecars).map(_ -> "--reference sidecar")
    protectedPaths.foreach { case (path, description) =>
      val samePath = output.toAbsolutePath.normalize == path.toAbsolutePath.normalize
      val sameFile = Files.exists(output) && Files.exists(path) && Files.isSameFile(output, path)
      require(!samePath && !sameFile, s"--output must not overwrite $description ($path), including file aliases")
    }
  }

  /** Write to a sibling and publish atomically; the writer must close successfully before returning. */
  private[tools] def writeOutput(write: Path => Unit): Unit = {
    val destination = output.toAbsolutePath
    // Io selects gzip/BGZF from the extension, so the temporary file must retain it.
    val temporary = Files.createTempFile(destination.getParent, ".rescue-junction-",
      PathUtil.extensionOf(output).getOrElse(".tmp"))
    try {
      write(temporary)
      assertOutputSafe()
      Files.move(temporary, destination, ATOMIC_MOVE, REPLACE_EXISTING)
    }
    finally Files.deleteIfExists(temporary)
  }

  override def execute(): Unit = {
    assertOutputSafe()
    val source = SamSource(input)
    try {
      require(source.header.getSortOrder == SortOrder.coordinate,
        "BAM must be coordinate-sorted; coordinate-sort and index --input")
      require(source.indexed, "BAM index is missing; index the coordinate-sorted --input BAM")
      val dictionary = source.dict
      val left = dictionary.get(leftContig).getOrElse(
        throw new IllegalArgumentException(s"Unknown --left-contig '$leftContig' in BAM; check candidate contigs"))
      val right = dictionary.get(rightContig).getOrElse(
        throw new IllegalArgumentException(s"Unknown --right-contig '$rightContig' in BAM; check candidate contigs"))
      require(leftPosition >= 1 && leftPosition <= left.length,
        s"--left-position must be in [1,${left.length}] for $leftContig")
      require(rightPosition >= 1 && rightPosition <= right.length,
        s"--right-position must be in [1,${right.length}] for $rightContig")
      require(ReferenceSequenceFileFactory.canCreateIndexedFastaReader(reference),
        "Reference FASTA index is missing or unusable; index --reference (and supply a .gzi for BGZF FASTA)")
      val fasta = ReferenceSequenceFileFactory.getReferenceSequenceFile(reference)
      try {
        require(fasta.isIndexed, "Reference FASTA must be indexed; index --reference")
        // The FASTA index carries lengths even when no sequence dictionary was supplied.
        val index = new FastaSequenceIndex(ReferenceSequenceFileFactory.getFastaIndexFileName(reference))
        IndexedSeq(left, right).distinct.foreach { contig =>
          require(index.hasIndexEntry(contig.name),
            s"Reference FASTA is missing BAM contig '${contig.name}'; use the BAM's matching --reference")
          require(index.getIndexEntry(contig.name).getSize == contig.length.toLong,
            s"Reference/BAM length mismatch for ${contig.name}; use the BAM's matching --reference")
          Option(fasta.getSequenceDictionary).foreach { referenceDictionary =>
            val entry = referenceDictionary.getSequence(contig.name)
            require(entry != null && entry.getSequenceLength == contig.length,
              s"Reference/BAM dictionary mismatch for ${contig.name}; use the BAM's matching --reference")
          }
        }
        val breakpoint = Breakpoint(left.index, leftPosition, leftStrand == "+",
          right.index, rightPosition, rightStrand == "+")
        val rescuer = new JunctionRescuer(
          breakpoint=breakpoint,
          reference=(index, start, end) => fasta.getSubsequenceAt(dictionary(index).name, start.toLong, end.toLong).getBases,
          contigLengths=dictionary.iterator.map(_.length).toIndexedSeq,
          flank=flank, minAnchor=minAnchor, slop=slop, minMapq=minMapq,
          minScoreDelta=minScoreDelta, minIdentity=minIdentity, minQueryCoverage=minQueryCoverage
        )
        val regions = IndexedSeq((left, leftPosition), (right, rightPosition)).map { case (contig, position) =>
          new Interval(contig.name, math.max(1L, position.toLong - flank).toInt,
            math.min(contig.length.toLong, position.toLong + flank).toInt)
        }
        val rescued = mutable.Map.empty[String, RescuedJunctionTemplate]
        // One multi-region query merges overlapping windows and visits each BAM record once.
        val iterator = source.query(regions, QueryType.Overlapping)
        try iterator.foreach { record =>
          rescuer.rescue(record).foreach { result =>
            val metric = RescuedJunctionTemplate(
              left_contig=leftContig, left_position=leftPosition, left_strand=leftStrand,
              right_contig=rightContig, right_position=rightPosition, right_strand=rightStrand,
              junction_boundary=rescuer.junctionBoundary,
              template=record.name, read_number=if (!record.paired) 0 else if (record.firstOfPair) 1 else 2,
              contig=record.refName, position=record.start, original_cigar=record.cigar.toString,
              score=result.score, reference_score=result.referenceScore,
              left_anchor=result.leftAnchor, right_anchor=result.rightAnchor,
              query_start=result.queryStart, query_end=result.queryEnd,
              target_start=result.targetStart, target_end=result.targetEnd, junction_cigar=result.cigar,
              query_reversed=result.reverseComplemented, identity=result.identity, query_coverage=result.queryCoverage
            )
            if (rescued.get(record.name).forall(previous => RescuedJunctionTemplate.representativeOrdering.lt(metric, previous))) {
              rescued(record.name) = metric
            }
          }
        }
        finally iterator.close()
        writeOutput { temporary =>
          Using.resource(Io.toWriter(temporary)) { writer =>
            Metric.write[RescuedJunctionTemplate](writer, rescued.values.toIndexedSeq.sortBy(_.template))
          }
        }
      }
      finally fasta.close()
    }
    finally source.close()
  }
}

/** Deterministic ordering within one template and candidate, with the best representative first. */
object RescuedJunctionTemplate {
  private[tools] val representativeOrdering: Ordering[RescuedJunctionTemplate] = Ordering.by { metric =>
    (
      (-(metric.score.toLong - metric.reference_score.toLong), -metric.score.toLong,
        -math.min(metric.left_anchor, metric.right_anchor).toLong),
      (metric.read_number, metric.contig, metric.position, metric.query_start, metric.query_end,
        metric.target_start, metric.target_end, (metric.original_cigar, metric.junction_cigar, metric.query_reversed)),
      (metric.left_anchor, metric.right_anchor, metric.identity, metric.query_coverage)
    )
  }
}

/** One qualifying template for the supplied junction; all coordinates are 1-based inclusive.
  *
  * @param left_contig supplied left/from breakend contig
  * @param left_position supplied left/from breakend position on the reference
  * @param left_strand supplied left/from breakend strand (+ or -)
  * @param right_contig supplied right/into breakend contig
  * @param right_position supplied right/into breakend position on the reference
  * @param right_strand supplied right/into breakend strand (+ or -)
  * @param junction_boundary actual incoming flank length after clamping; the last left-side target base
  * @param template BAM query/template name; existing discovery supporters are not excluded
  * @param read_number representative read number: 0 for unpaired, 1 for first, 2 for second
  * @param contig representative original BAM alignment contig
  * @param position representative original BAM alignment start
  * @param original_cigar representative original BAM CIGAR
  * @param score junction Smith-Waterman alignment score
  * @param reference_score best score against either local reference alternative in either query orientation
  * @param left_anchor canonical aligned columns on the incoming side, including mismatches but not gaps/N
  * @param right_anchor canonical aligned columns on the outgoing side, including mismatches but not gaps/N
  * @param query_start alignment start on the query orientation indicated by query_reversed
  * @param query_end alignment end on the query orientation indicated by query_reversed
  * @param target_start alignment start on the oriented synthetic junction
  * @param target_end alignment end on the oriented synthetic junction
  * @param junction_cigar realigned CIGAR against the synthetic junction, not a replacement BAM CIGAR
  * @param query_reversed true if query coordinates use reverse-complemented SAM-stored bases
  * @param identity fraction of aligned columns that are canonical matching bases
  * @param query_coverage fraction of observed query bases included in the junction alignment
  */
case class RescuedJunctionTemplate(left_contig: String,
                                  left_position: Int,
                                  left_strand: String,
                                  right_contig: String,
                                  right_position: Int,
                                  right_strand: String,
                                  junction_boundary: Int,
                                  template: String,
                                  read_number: Int,
                                  contig: String,
                                  position: Int,
                                  original_cigar: String,
                                  score: Int,
                                  reference_score: Int,
                                  left_anchor: Int,
                                  right_anchor: Int,
                                  query_start: Int,
                                  query_end: Int,
                                  target_start: Int,
                                  target_end: Int,
                                  junction_cigar: String,
                                  query_reversed: Boolean,
                                  identity: Double,
                                  query_coverage: Double) extends Metric
