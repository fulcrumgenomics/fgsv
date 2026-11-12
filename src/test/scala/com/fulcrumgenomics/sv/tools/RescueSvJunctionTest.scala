package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.bam.api.{SamOrder, SamWriter}
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.sv.UnitSpec
import com.fulcrumgenomics.sv.cmdline.SvMain
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder}
import com.fulcrumgenomics.util.{Io, Metric}
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

import java.io.{FilterWriter, IOException}
import java.nio.file.{Files, Path}
import scala.jdk.CollectionConverters._
import scala.util.Using

class RescueSvJunctionTest extends UnitSpec {
  private val sequence = {
    val random = new scala.util.Random(42)
    (1 to 1000).map(_ => "ACGT"(random.nextInt(4))).mkString
  }
  private val bases = sequence.substring(575, 600) + sequence.substring(400, 525)

  private def reference(name: String = "chr1", length: Int = 1000, dictionary: Boolean = true): Path = {
    val refs = new ReferenceSetBuilder
    refs.add(name).add(sequence.take(length))
    val fasta = refs.toTempFile()
    if (!dictionary) Files.delete(ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta))
    fasta
  }

  private def builder(readLength: Int = 150, sort: SamOrder = SamOrder.Coordinate): SamBuilder = {
    new SamBuilder(readLength=readLength, sort=Some(sort),
      sd=Some(SequenceDictionary(SequenceMetadata("chr1", 1000))))
  }

  private def tool(input: Path,
                   fasta: Path,
                   output: Path,
                   leftContig: String = "chr1",
                   leftPosition: Int = 600,
                   leftStrand: String = "+",
                   rightContig: String = "chr1",
                   rightPosition: Int = 401,
                   rightStrand: String = "+",
                   minAnchor: Int = 10,
                   slop: Int = 5,
                   minMapq: Int = 30,
                   minScoreDelta: Int = 10,
                   minIdentity: Double = 0.9,
                   minQueryCoverage: Double = 0.9,
                   flank: Int = 250): RescueSvJunction = {
    new RescueSvJunction(input=input, reference=fasta, output=output,
      leftContig=leftContig, leftPosition=leftPosition, leftStrand=leftStrand,
      rightContig=rightContig, rightPosition=rightPosition, rightStrand=rightStrand,
      minAnchor=minAnchor, slop=slop, minMapq=minMapq, minScoreDelta=minScoreDelta,
      minIdentity=minIdentity, minQueryCoverage=minQueryCoverage, flank=flank)
  }

  /** Enumerate protected sidecars independently of the production output-safety check. */
  private def sidecars(bam: Path, fasta: Path): IndexedSeq[(Path, String)] = {
    val basename = bam.resolveSibling(bam.getFileName.toString.stripSuffix(".bam"))
    val bamPaths = IndexedSeq(".bai", ".csi").flatMap { suffix =>
      IndexedSeq(bam, basename).map(path => path.resolveSibling(path.getFileName.toString + suffix))
    }
    val referencePaths = IndexedSeq(
      ReferenceSequenceFileFactory.getFastaIndexFileName(fasta),
      fasta.resolveSibling(fasta.getFileName.toString + ".gzi"),
      ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta)
    )
    bamPaths.map(_ -> "--input index") ++ referencePaths.map(_ -> "--reference sidecar")
  }

  /** Isolate publication tests from unrelated temporary files created by concurrent suites. */
  private def isolatedOutput(suffix: String): Path = {
    val directory = Io.makeTempDir("rescue-output-")
    directory.toFile.deleteOnExit()
    val output = Files.createTempFile(directory, "metrics-", suffix)
    output.toFile.deleteOnExit()
    output
  }

  /** List directory entries so tests can detect leaked temporary siblings. */
  private def entries(directory: Path): Set[Path] = {
    Using.resource(Files.list(directory))(_.iterator().asScala.toSet)
  }

  /** Run a failing command against an absent output, checking that it did not emit a metrics header. */
  private def failsWithoutOutput(fragment: String)(run: Path => Unit): Unit = {
    val output = makeTempFile("failed-rescue", ".txt")
    Files.delete(output)
    val error = intercept[IllegalArgumentException] { run(output) }
    error.getMessage should include (fragment)
    Files.exists(output) shouldBe false
  }

  "RescueSvJunction" should "rescue real paired reads once per template without changing any input bytes" in {
    val fasta = reference()
    val reads = builder()
    // A tandem duplication joins chr1:600 to chr1:401; the best mate is selected, not counted twice.
    reads.addPair(name="rescued", start1=401, start2=401, cigar1="25S125M", cigar2="10S140M",
      bases1=bases, bases2=sequence.substring(590, 600) + sequence.substring(400, 540),
      strand1=SamBuilder.Plus, strand2=SamBuilder.Plus)
    reads.addFrag(name="ten-base", start=401, cigar="10S140M",
      bases=sequence.substring(590, 600) + sequence.substring(400, 540))
    reads.addFrag(name="nine-base", start=401, cigar="9S141M",
      bases=sequence.substring(591, 600) + sequence.substring(400, 541))
    reads.addFrag(name="short", start=401, cigar="3S147M",
      bases=sequence.substring(597, 600) + sequence.substring(400, 547))
    reads.addFrag(name="reference", start=401, cigar="25S125M", bases=sequence.substring(375, 525))
    val bam = reads.toTempFile()
    val original = Files.readAllBytes(bam)
    val originalFasta = Files.readAllBytes(fasta)
    val output = makeTempFile("rescued", ".txt")
    tool(bam, fasta, output).execute()
    val rows = Metric.read[RescuedJunctionTemplate](output).toIndexedSeq
    rows.map(_.template) shouldBe IndexedSeq("rescued", "ten-base")
    val row = rows.head
    row.read_number shouldBe 1
    row.left_anchor shouldBe 25
    row.right_anchor shouldBe 125
    row.score shouldBe 300
    row.reference_score shouldBe 250
    (row.left_contig, row.left_position, row.left_strand) shouldBe (("chr1", 600, "+"))
    (row.right_contig, row.right_position, row.right_strand) shouldBe (("chr1", 401, "+"))
    row.junction_boundary shouldBe 250
    row.target_start shouldBe row.junction_boundary - 24
    row.target_end shouldBe row.junction_boundary + 125
    (row.query_start, row.query_end, row.query_reversed) shouldBe ((1, 150, false))
    Files.readAllBytes(bam).toIndexedSeq shouldBe original.toIndexedSeq
    Files.readAllBytes(fasta).toIndexedSeq shouldBe originalFasta.toIndexedSeq
    val pair = readBamRecs(bam).filter(_.name == "rescued")
    pair.size shouldBe 2
    pair.forall(_.paired) shouldBe true
    pair.map(_.firstOfPair).toSet shouldBe Set(true, false)
    pair.foreach(_.get[String]("SA") shouldBe None)
  }

  it should "produce identical rescue metrics for literal and SAM equals bases in indexed BAMs" in {
    val fasta = reference()
    val outputs = IndexedSeq(bases, bases.take(25) + ("=" * 125)).map { storedBases =>
      val reads = builder()
      reads.addFrag(name="rescued", start=401, cigar="25S125M", bases=storedBases)
      val bam = reads.toTempFile()
      readBamRecs(bam).head.basesString shouldBe storedBases
      val original = Files.readAllBytes(bam).toIndexedSeq
      val output = isolatedOutput(".txt")
      tool(bam, fasta, output).execute()
      val metrics = Metric.read[RescuedJunctionTemplate](output).toIndexedSeq
      metrics.map(_.template) shouldBe IndexedSeq("rescued")
      Files.readAllBytes(bam).toIndexedSeq shouldBe original
      Files.readAllBytes(output).toIndexedSeq
    }
    outputs.head shouldBe outputs.last
  }

  it should "report the actual clamped junction boundary without requiring a FASTA dictionary" in {
    val fasta = reference(dictionary=false)
    val reads = builder()
    reads.addFrag(name="edge", start=401, cigar="25S125M",
      bases=sequence.substring(0, 25) + sequence.substring(400, 525))
    val output = makeTempFile("edge-rescue", ".txt")
    tool(reads.toTempFile(), fasta, output, leftPosition=25, flank=200).execute()
    val rows = Metric.read[RescuedJunctionTemplate](output).toIndexedSeq
    rows.size shouldBe 1
    val row = rows.head
    row.junction_boundary shouldBe 25
    row.target_start shouldBe 1
    row.target_end shouldBe 150
    row.left_position shouldBe 25
  }

  it should "produce the same representative and bytes for tied paired reads in either input order" in {
    val fasta = reference()
    val reads = builder()
    val pair = reads.addPair(name="tied", start1=401, start2=401, cigar1="25S125M", cigar2="25S125M",
      bases1=bases, bases2=bases, strand1=SamBuilder.Plus, strand2=SamBuilder.Plus)
    // Both query windows overlap this pair. The multi-region query must still emit one row.
    val outputs = IndexedSeq(pair, pair.reverse).map { records =>
      val bam = makeTempFile("input-order", ".bam")
      val writer = SamWriter(bam, reads.header)
      try writer ++= records finally writer.close()
      val output = makeTempFile("tied-rescue", ".txt")
      tool(bam, fasta, output, flank=200).execute()
      val rows = Metric.read[RescuedJunctionTemplate](output).toIndexedSeq
      rows.size shouldBe 1
      rows.head.read_number shouldBe 1
      Files.readAllBytes(output).toIndexedSeq
    }
    outputs.head shouldBe outputs.last
  }

  it should "totally order varying diagnostics within one template and candidate by the documented ranking" in {
    val best = RescuedJunctionTemplate("chr1", 600, "+", "chr1", 401, "+", 200,
      "template", 1, "chr1", 401, "25S125M", 300, 250, 25, 125, 1, 150, 176, 325, "150M", false, 1.0, 1.0)
    val ordering = RescuedJunctionTemplate.representativeOrdering
    val ranked = IndexedSeq(
      best,
      best.copy(read_number=2),
      best.copy(read_number=0, left_anchor=24),
      best.copy(read_number=0, score=299, reference_score=249),
      best.copy(read_number=0, score=310, reference_score=261)
    )
    ranked.reverse.sorted(ordering) shouldBe ranked
    val ascendingTieBreaks = IndexedSeq(
      best.copy(read_number=2), best.copy(contig="chr2"), best.copy(position=402),
      best.copy(query_start=2), best.copy(query_end=151), best.copy(target_start=177), best.copy(target_end=326),
      best.copy(original_cigar="26S124M"), best.copy(junction_cigar="151M"), best.copy(query_reversed=true)
    )
    ascendingTieBreaks.foreach { other =>
      ordering.lt(best, other) shouldBe true
      IndexedSeq(other, best).sorted(ordering) shouldBe IndexedSeq(best, other)
    }
    val diagnosticTieBreaks = IndexedSeq(best.copy(right_anchor=126), best.copy(left_anchor=26),
      best.copy(identity=0.99), best.copy(query_coverage=0.99))
    (ascendingTieBreaks ++ diagnosticTieBreaks).foreach { other =>
      ordering.compare(best, other) should not be 0
      ordering.compare(best, other) shouldBe -ordering.compare(other, best)
      IndexedSeq(other, best).sorted(ordering) shouldBe IndexedSeq(best, other).sorted(ordering)
    }
  }

  it should "write the same complete header for empty input and input with no qualifying templates" in {
    val fasta = reference()
    val empty = builder().toTempFile()
    val reads = builder()
    reads.addFrag(name="reference", start=401, cigar="25S125M", bases=sequence.substring(375, 525))
    val headers = IndexedSeq(empty, reads.toTempFile()).map { bam =>
      val output = makeTempFile("empty-rescue", ".txt")
      tool(bam, fasta, output).execute()
      Metric.read[RescuedJunctionTemplate](output).toIndexedSeq shouldBe IndexedSeq.empty
      val lines = Io.readLines(output).toIndexedSeq
      lines.size shouldBe 1
      lines.head.split('\t').toIndexedSeq
    }
    headers.head shouldBe headers.last
    headers.head shouldBe IndexedSeq("left_contig", "left_position", "left_strand", "right_contig", "right_position",
      "right_strand", "junction_boundary", "template", "read_number", "contig", "position", "original_cigar", "score",
      "reference_score", "left_anchor", "right_anchor", "query_start", "query_end", "target_start", "target_end",
      "junction_cigar", "query_reversed", "identity", "query_coverage")
  }

  it should "preserve existing output or leave no output and clean siblings when a writer fails after partial output" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    IndexedSeq(".txt", ".gz", ".bgz", ".bgzip").foreach { suffix =>
      IndexedSeq(false, true).foreach { existing =>
        val output = isolatedOutput(suffix)
        val sentinel = IndexedSeq[Byte](11, 22, 33)
        if (existing) Files.write(output, sentinel.toArray) else Files.delete(output)
        val command = tool(bam, fasta, output)
        val before = entries(output.getParent)
        val failure = new IOException("injected late write failure")
        val error = intercept[IOException] {
          command.writeOutput { temporary =>
            temporary.getParent shouldBe output.getParent
            temporary should not be output
            temporary.getFileName.toString should endWith (suffix)
            Using.resource(Io.toWriter(temporary)) { writer =>
              writer.write("partial metrics output")
              writer.flush()
              Files.size(temporary) should be > 0L
              throw failure
            }
          }
        }
        error shouldBe failure
        Files.exists(output) shouldBe existing
        if (existing) Files.readAllBytes(output).toIndexedSeq shouldBe sentinel
        entries(output.getParent) shouldBe before
      }
    }
  }

  it should "propagate a Metric writer close failure without publishing output or leaving temporary siblings" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    IndexedSeq(".txt", ".gz", ".bgz", ".bgzip").foreach { suffix =>
      IndexedSeq(false, true).foreach { existing =>
        val output = isolatedOutput(suffix)
        val sentinel = IndexedSeq[Byte](11, 22, 33)
        if (existing) Files.write(output, sentinel.toArray) else Files.delete(output)
        val command = tool(bam, fasta, output)
        val before = entries(output.getParent)
        val failure = new IOException("injected close failure")
        var closeAttempted = false
        val error = intercept[IOException] {
          command.writeOutput { temporary =>
            val failingWriter = new FilterWriter(Io.toWriter(temporary)) {
              override def close(): Unit = {
                super.close()
                if (!closeAttempted) {
                  closeAttempted = true
                  // Fail only on the first close, so a second close cannot mask a swallowed failure.
                  Files.size(temporary) should be > 0L
                  throw failure
                }
              }
            }
            // Match execute's resource handling and Metric.write overload exactly.
            Using.resource(failingWriter) { writer =>
              Metric.write[RescuedJunctionTemplate](writer, IndexedSeq.empty[RescuedJunctionTemplate])
            }
          }
        }
        error shouldBe failure
        closeAttempted shouldBe true
        Files.exists(output) shouldBe existing
        if (existing) Files.readAllBytes(output).toIndexedSeq shouldBe sentinel
        entries(output.getParent) shouldBe before
      }
    }
  }

  it should "publish readable plain, gzip, and BGZF metrics without leaving temporary siblings" in {
    val fasta = reference()
    val reads = builder()
    reads.addFrag(name="rescued", start=401, cigar="25S125M", bases=bases)
    val bam = reads.toTempFile()
    val rows = IndexedSeq(".txt", ".gz", ".bgz", ".bgzip").flatMap { suffix =>
      IndexedSeq(false, true).map { existing =>
        val output = isolatedOutput(suffix)
        if (existing) Files.write(output, Array[Byte](11, 22, 33)) else Files.delete(output)
        val before = entries(output.getParent)
        tool(bam, fasta, output).execute()
        entries(output.getParent) shouldBe (before + output)
        val metrics = Metric.read[RescuedJunctionTemplate](output).toIndexedSeq
        metrics.map(_.template) shouldBe IndexedSeq("rescued")
        val bytes = Files.readAllBytes(output)
        if (suffix != ".txt") {
          bytes.take(2).toIndexedSeq shouldBe IndexedSeq[Byte](31, -117)
          // BGZF uses the gzip extra-field flag, unlike Io's ordinary gzip output.
          bytes(3) shouldBe (if (suffix == ".gz") 0.toByte else 4.toByte)
        }
        metrics
      }
    }
    rows.distinct.size shouldBe 1
  }

  it should "recheck input and sidecar aliases immediately before publication and clean temporary siblings" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    val protectedPaths = IndexedSeq(bam -> "--input", fasta -> "--reference") ++
      sidecars(bam, fasta).filter { case (path, _) => Files.exists(path) }
    protectedPaths.foreach { case (path, description) =>
      IndexedSeq(false, true).foreach { symlink =>
        val output = isolatedOutput(".txt")
        val command = tool(bam, fasta, output)
        val original = Files.readAllBytes(path).toIndexedSeq
        val before = entries(output.getParent)
        val error = intercept[IllegalArgumentException] {
          command.writeOutput { temporary =>
            Io.writeLines(temporary, IndexedSeq("complete temporary output"))
            Files.delete(output)
            if (symlink) Files.createSymbolicLink(output, path) else Files.createLink(output, path)
          }
        }
        error.getMessage should include (description)
        Files.isSameFile(output, path) shouldBe true
        Files.readAllBytes(path).toIndexedSeq shouldBe original
        entries(output.getParent) shouldBe before
      }
    }
  }

  it should "clean temporary siblings when atomic publication fails" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    val output = isolatedOutput(".txt")
    val command = tool(bam, fasta, output)
    Files.delete(output)
    Files.createDirectory(output)
    val child = output.resolve("sentinel")
    Files.write(child, Array[Byte](11, 22, 33))
    try {
      val before = entries(output.getParent)
      intercept[IOException] {
        command.writeOutput(temporary => Io.writeLines(temporary, IndexedSeq("complete temporary output")))
      }
      Files.readAllBytes(child).toIndexedSeq shouldBe IndexedSeq[Byte](11, 22, 33)
      entries(output.getParent) shouldBe before
    }
    finally {
      Files.delete(child)
      Files.delete(output)
    }
  }

  it should "parse command-line defaults and equals-form negative strands" in {
    val fasta = reference()
    val bam = builder().toTempFile()
    val output = makeTempFile("parsed-rescue", ".txt")
    new SvMain().makeItSo(Array("RescueSvJunction", "--input", bam.toString,
      "--reference", fasta.toString, "--output", output.toString,
      "--left-contig", "chr1", "--left-position", "600", "--left-strand=-",
      "--right-contig", "chr1", "--right-position", "401", "--right-strand=-")) shouldBe 0
    Io.readLines(output).size shouldBe 1
    Metric.read[RescuedJunctionTemplate](output).toIndexedSeq shouldBe IndexedSeq.empty
  }

  it should "ignore oversize ineligible records but reject an otherwise eligible read that exceeds flank minus slop" in {
    val fasta = reference()
    val reads = builder(readLength=200)
    val longBases = sequence.substring(575, 600) + sequence.substring(400, 575)
    reads.addFrag(name="duplicate", start=401, cigar="25S175M", bases=longBases).foreach(_.duplicate = true)
    reads.addFrag(name="secondary", start=401, cigar="25S175M", bases=longBases).foreach(_.secondary = true)
    reads.addFrag(name="supplementary", start=401, cigar="25S175M", bases=longBases).foreach(_.supplementary = true)
    reads.addFrag(name="failed-qc", start=401, cigar="25S175M", bases=longBases).foreach(_.pf = false)
    reads.addFrag(name="low-mapq", start=401, cigar="25S175M", bases=longBases, mapq=0)
    reads.addFrag(name="unclipped", start=401, cigar="200M", bases=sequence.substring(400, 600))
    reads.addFrag(name="distant-clip", start=450, cigar="25S175M", bases=longBases)
    val bam = reads.toTempFile()
    val original = Files.readAllBytes(bam).toIndexedSeq
    val output = makeTempFile("ineligible", ".txt")
    tool(bam, fasta, output, flank=200).execute()
    Metric.read[RescuedJunctionTemplate](output).toIndexedSeq shouldBe IndexedSeq.empty
    Files.readAllBytes(bam).toIndexedSeq shouldBe original

    reads.addFrag(name="eligible-long", start=401, cigar="25S175M", bases=longBases)
    val eligibleBam = reads.toTempFile()
    val before = Files.readAllBytes(eligibleBam).toIndexedSeq
    failsWithoutOutput("flank") { failedOutput => tool(eligibleBam, fasta, failedOutput, flank=200).execute() }
    Files.readAllBytes(eligibleBam).toIndexedSeq shouldBe before
    val sufficientOutput = makeTempFile("sufficient-flank", ".txt")
    tool(eligibleBam, fasta, sufficientOutput, flank=205).execute()
    Metric.read[RescuedJunctionTemplate](sufficientOutput).toIndexedSeq.map(_.template) shouldBe IndexedSeq("eligible-long")
  }

  it should "leave output untouched when an eligible oversize read fails after a qualifying template was collected" in {
    val fasta = reference()
    val reads = builder()
    reads.addFrag(name="qualifying-first", start=401, cigar="25S125M", bases=bases)
    val qualifyingBam = reads.toTempFile()
    val successfulOutput = isolatedOutput(".txt")
    tool(qualifyingBam, fasta, successfulOutput, flank=200).execute()
    Metric.read[RescuedJunctionTemplate](successfulOutput).toIndexedSeq.map(_.template) shouldBe IndexedSeq("qualifying-first")

    // A distinct, later coordinate guarantees the qualifying read is encountered before this failure.
    reads.addFrag(name="oversize-later", start=402, cigar="25S175M",
      bases=sequence.substring(575, 600) + sequence.substring(401, 576), quals="I" * 200)
    val bam = reads.toTempFile()
    readBamRecs(bam).map(_.name) shouldBe IndexedSeq("qualifying-first", "oversize-later")
    val original = Files.readAllBytes(bam).toIndexedSeq
    IndexedSeq(false, true).foreach { existing =>
      val output = isolatedOutput(".txt")
      val sentinel = IndexedSeq[Byte](11, 22, 33)
      if (existing) Files.write(output, sentinel.toArray) else Files.delete(output)
      val before = entries(output.getParent)
      intercept[IllegalArgumentException] { tool(bam, fasta, output, flank=200).execute() }
        .getMessage should include ("flank")
      Files.exists(output) shouldBe existing
      if (existing) Files.readAllBytes(output).toIndexedSeq shouldBe sentinel
      entries(output.getParent) shouldBe before
    }
    Files.readAllBytes(bam).toIndexedSeq shouldBe original
  }

  it should "reject invalid strands and thresholds with actionable errors before writing output" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    val invalid = IndexedSeq[(String, Path => RescueSvJunction)](
      "--left-strand" -> (out => tool(bam, fasta, out, leftStrand="?")),
      "--right-strand" -> (out => tool(bam, fasta, out, rightStrand="?")),
      "--min-anchor" -> (out => tool(bam, fasta, out, minAnchor=0)),
      "--slop" -> (out => tool(bam, fasta, out, slop= -1)),
      "--flank" -> (out => tool(bam, fasta, out, flank=5)),
      "--flank" -> (out => tool(bam, fasta, out, flank=9, minAnchor=10, slop=0)),
      "--min-mapq" -> (out => tool(bam, fasta, out, minMapq= -1)),
      "--min-mapq" -> (out => tool(bam, fasta, out, minMapq=256)),
      "--min-score-delta" -> (out => tool(bam, fasta, out, minScoreDelta= -1)),
      "--min-identity" -> (out => tool(bam, fasta, out, minIdentity=0)),
      "--min-identity" -> (out => tool(bam, fasta, out, minIdentity=1.1)),
      "--min-identity" -> (out => tool(bam, fasta, out, minIdentity=Double.NaN)),
      "--min-query-coverage" -> (out => tool(bam, fasta, out, minQueryCoverage=0)),
      "--min-query-coverage" -> (out => tool(bam, fasta, out, minQueryCoverage=1.1)),
      "--min-query-coverage" -> (out => tool(bam, fasta, out, minQueryCoverage=Double.NaN))
    )
    invalid.foreach { case (fragment, makeTool) => failsWithoutOutput(fragment)(out => { makeTool(out); () }) }
  }

  it should "reject unknown candidate contigs and out-of-range coordinates before writing output" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    failsWithoutOutput("--left-contig") { out => tool(bam, fasta, out, leftContig="missing").execute() }
    failsWithoutOutput("--right-contig") { out => tool(bam, fasta, out, rightContig="missing").execute() }
    IndexedSeq(0, 1001).foreach { position =>
      failsWithoutOutput("--left-position") { out => tool(bam, fasta, out, leftPosition=position).execute() }
      failsWithoutOutput("--right-position") { out => tool(bam, fasta, out, rightPosition=position).execute() }
    }
  }

  it should "reject missing BAM and FASTA indices and non-coordinate BAMs before writing output" in {
    val fasta = reference()
    val reads = builder()
    val bam = makeTempFile("unindexed", ".bam")
    val writer = SamWriter(bam, reads.header, index=false)
    writer.close()
    failsWithoutOutput("BAM index") { out => tool(bam, fasta, out).execute() }
    val unsorted = builder(sort=SamOrder.Unsorted).toTempFile()
    failsWithoutOutput("coordinate-sort") { out => tool(unsorted, fasta, out).execute() }
    val indexedBam = builder().toTempFile()
    Files.delete(ReferenceSequenceFileFactory.getFastaIndexFileName(fasta))
    failsWithoutOutput("FASTA index") { out => tool(indexedBam, fasta, out).execute() }
  }

  it should "reject reference contig and length mismatches even without a FASTA dictionary" in {
    val bam = builder().toTempFile()
    IndexedSeq(true, false).foreach { dictionary =>
      val wrongContig = reference(name="other", dictionary=dictionary)
      failsWithoutOutput("missing BAM contig 'chr1'") { out => tool(bam, wrongContig, out).execute() }
      val wrongLength = reference(length=999, dictionary=dictionary)
      failsWithoutOutput("length mismatch for chr1") { out => tool(bam, wrongLength, out).execute() }
    }
  }

  it should "refuse to overwrite inputs or sidecars through direct paths, symlinks, or hardlinks" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    val protectedPaths = IndexedSeq(bam -> "--input", fasta -> "--reference") ++ sidecars(bam, fasta)
    protectedPaths.foreach { case (path, _) =>
      if (!Files.exists(path)) {
        Files.write(path, Array[Byte](1, 2, 3))
        path.toFile.deleteOnExit()
      }
    }
    val before = protectedPaths.map { case (path, _) => path -> Files.readAllBytes(path).toIndexedSeq }
    protectedPaths.foreach { case (path, description) =>
      val symlink = makeTempFile("symlink-output", ".txt")
      Files.delete(symlink)
      Files.createSymbolicLink(symlink, path)
      val hardlink = makeTempFile("hardlink-output", ".txt")
      Files.delete(hardlink)
      Files.createLink(hardlink, path)
      IndexedSeq(path, symlink, hardlink).foreach { output =>
        intercept[IllegalArgumentException] { tool(bam, fasta, output).execute() }.getMessage should include (description)
        Files.readAllBytes(output).toIndexedSeq shouldBe Files.readAllBytes(path).toIndexedSeq
      }
    }
    before.foreach { case (path, bytes) => Files.readAllBytes(path).toIndexedSeq shouldBe bytes }
  }

  it should "reserve missing input sidecar paths rather than write metrics into them" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    val beforeBam = Files.readAllBytes(bam).toIndexedSeq
    val beforeFasta = Files.readAllBytes(fasta).toIndexedSeq
    sidecars(bam, fasta).foreach { case (output, _) =>
      Files.deleteIfExists(output)
      intercept[IllegalArgumentException] { tool(bam, fasta, output) }.getMessage should include ("--output must not overwrite")
      Files.exists(output) shouldBe false
    }
    Files.readAllBytes(bam).toIndexedSeq shouldBe beforeBam
    Files.readAllBytes(fasta).toIndexedSeq shouldBe beforeFasta
  }

  it should "reject dangling output symlinks to absent dictionaries and BAM sidecars before execution or publication" in {
    val bam = builder().toTempFile()
    val fasta = reference(dictionary=false)
    // Keep the usable BAM index; reserve the other conventional index spellings too.
    val missing = sidecars(bam, fasta).map(_._1).filterNot(Files.exists(_))
    missing should contain (ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(fasta))
    missing.exists(_.getFileName.toString.endsWith(".csi")) shouldBe true
    val protectedPaths = IndexedSeq(bam, fasta) ++ sidecars(bam, fasta).map(_._1).filter(Files.exists(_))
    val original = protectedPaths.map(path => path -> Files.readAllBytes(path).toIndexedSeq)
    missing.foreach { target =>
      IndexedSeq(false, true).foreach { relative =>
        IndexedSeq("construction", "execution", "publication").foreach { phase =>
          val output = isolatedOutput(".txt")
          val linkTarget = if (relative) output.getParent.relativize(target) else target.toAbsolutePath
          def introduceLink(): Unit = {
            Files.delete(output)
            Files.createSymbolicLink(output, linkTarget)
          }
          val command = if (phase == "construction") {
            introduceLink()
            None
          }
          else Some(tool(bam, fasta, output))
          if (phase == "execution") introduceLink()
          val before = entries(output.getParent)
          val error = intercept[IllegalArgumentException] {
            phase match {
              case "construction" => tool(bam, fasta, output)
              case "execution" => command.value.execute()
              case "publication" => command.value.writeOutput { temporary =>
                Io.writeLines(temporary, IndexedSeq("complete temporary output"))
                introduceLink()
              }
            }
          }
          error.getMessage should include ("dangling symbolic link")
          Files.isSymbolicLink(output) shouldBe true
          Files.readSymbolicLink(output) shouldBe linkTarget
          Files.exists(target) shouldBe false
          entries(output.getParent) shouldBe before
          original.foreach { case (path, bytes) => Files.readAllBytes(path).toIndexedSeq shouldBe bytes }
        }
      }
    }
  }

  it should "recheck output identity before execution when an output alias is introduced after construction" in {
    val bam = builder().toTempFile()
    val fasta = reference()
    val before = Files.readAllBytes(bam).toIndexedSeq
    val output = makeTempFile("late-alias", ".txt")
    val command = tool(bam, fasta, output)
    Files.delete(output)
    Files.createLink(output, bam)
    intercept[IllegalArgumentException] { command.execute() }.getMessage should include ("--input")
    Files.readAllBytes(bam).toIndexedSeq shouldBe before
  }
}
