package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.FgBioDef.{FilePath, PathPrefix, PathToBam}
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.{SamRecord, SamWriter}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.sv.EvidenceType.{ReadPair, SplitRead}
import com.fulcrumgenomics.sv.SegmentOrigin.{Both, ReadOne, ReadTwo}
import com.fulcrumgenomics.sv._
import com.fulcrumgenomics.testing.SamBuilder
import htsjdk.samtools.SamPairUtil

class SvPileupTest extends UnitSpec {
  private def fromRangeOnly(refIndex: Int, start: Int, end: Int, origin: SegmentOrigin = ReadOne, positive: Boolean = true): AlignedSegment =
    AlignedSegment(
      origin=origin,
      readStart=1, readEnd=100, positiveStrand=positive, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=refIndex, start=start, end=end)
    )

  "SvPileup.isInterContigBreakpoint" should "is if two alignment segments are an inter-contig breakpoint" in {
    val chr1 = fromRangeOnly(refIndex=0, start=1, end=100)
    val chr2 = fromRangeOnly(refIndex=1, start=1, end=100)

    // Same contig
    SvPileup.isInterContigBreakpoint(seg1=chr1, seg2=chr1) shouldBe false

    // Different contig
    SvPileup.isInterContigBreakpoint(seg1=chr1, seg2=chr2) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr2, seg2=chr1) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr1, seg2=chr2.copy(origin=ReadTwo)) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr2.copy(origin=ReadTwo), seg2=chr1) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr1, seg2=chr2.copy(origin=Both)) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr2.copy(origin=Both), seg2=chr1) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr1.copy(origin=Both), seg2=chr2.copy(origin=Both)) shouldBe true
    SvPileup.isInterContigBreakpoint(seg1=chr2.copy(origin=Both), seg2=chr1.copy(origin=Both)) shouldBe true
  }

  "SvPileup.isIntraContigBreakpoint" should "identify where two segments are not moving up the genome within the defined spacing" in {
    val earlier = fromRangeOnly(refIndex=1, start=1, end=100)
    val overlap = fromRangeOnly(refIndex=1, start=50, end=150)
    val later   = fromRangeOnly(refIndex=1, start=150, end=200)

    // It has to be recalled that this test is done _after_ overlapping segments for R1 and R2 have been merged so
    // any overlapping segments or jumping backwards between segments is indicative of a breakpoint

    // What you might get from a read pair with a gap between the two reads
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=ReadTwo), maxWithinReadDistance=5, maxBetweenReadDistance=150) shouldBe false

    // same segment but jumping backwards
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=earlier, maxWithinReadDistance=5, maxBetweenReadDistance=150) shouldBe true

    // overlapping segments but still jumping backwards
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=overlap, maxWithinReadDistance=5, maxBetweenReadDistance=150) shouldBe true

    // non-overlapping segments from the same read, testing various values for maxWithinReadDistance
    // Note that the inner distance between two blocks is defined as `later.start - earlier.end`, so for
    // this case that is 150-100 = 50, so a breakpoint should be called when the maxWithinReadDistance < 50.
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later, maxWithinReadDistance= 5, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later, maxWithinReadDistance=25, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later, maxWithinReadDistance=48, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later, maxWithinReadDistance=49, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later, maxWithinReadDistance=50, maxBetweenReadDistance=150) shouldBe false
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later, maxWithinReadDistance=51, maxBetweenReadDistance=150) shouldBe false

    // non-overlapping segments where the later segment is "both" so indicates a split read breakpoint
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=Both), maxWithinReadDistance= 5, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=Both), maxWithinReadDistance=25, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=Both), maxWithinReadDistance=49, maxBetweenReadDistance=150) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=Both), maxWithinReadDistance=75, maxBetweenReadDistance=150) shouldBe false

    // non-overlapping segments where the later segment is "Read2" so between read distance should be used
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=ReadTwo), maxWithinReadDistance=5, maxBetweenReadDistance=5 ) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=ReadTwo), maxWithinReadDistance=5, maxBetweenReadDistance=25) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=ReadTwo), maxWithinReadDistance=5, maxBetweenReadDistance=49) shouldBe true
    SvPileup.isIntraContigBreakpoint(seg1=earlier, seg2=later.copy(origin=ReadTwo), maxWithinReadDistance=5, maxBetweenReadDistance=75) shouldBe false
  }

  "SvPileup.isIntraContigBreakpoint" should "identify when two segments flip strand" in {
    // By the time we're in segment space, all segments from a well-formed FR read pair should
    // share the same strand, so only when strands differ should a breakpoint be called

    val f1 = fromRangeOnly(refIndex=1, start=1,   end=99)
    val f2 = fromRangeOnly(refIndex=1, start=200, end=299)
    val r1 = f1.copy(positiveStrand=false)
    val r2 = f2.copy(positiveStrand=false)

    // Simple tests that should not call breakpoints
    SvPileup.isIntraContigBreakpoint(f1, f2, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe false
    SvPileup.isIntraContigBreakpoint(r2, r1, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe false

    // Now what if we make them different reads
    SvPileup.isIntraContigBreakpoint(f1, f2.copy(origin=ReadTwo), maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe false
    SvPileup.isIntraContigBreakpoint(r2, r1.copy(origin=ReadTwo), maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe false

    // But any combination on different strands should yield a breakpoint
    SvPileup.isIntraContigBreakpoint(f1, r1, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(f1, r2, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(f2, r1, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(f2, r2, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(r1, f1, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(r1, f2, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(r2, f1, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
    SvPileup.isIntraContigBreakpoint(r2, f2, maxWithinReadDistance=500, maxBetweenReadDistance=500) shouldBe true
  }

  //////////////////////////////////////////////////////////////////////////////
  // Objects and functions used in testing findBreakpoint()
  //////////////////////////////////////////////////////////////////////////////
  private val builder = new SamBuilder(readLength=100)
  import SamBuilder.{Minus, Plus, Strand}

  /** Construct a read/rec with the information necessary for breakpoint detection. */
  def r(chrom: String, pos: Int, strand: SamBuilder.Strand, r: Int, cigar: String, supp: Boolean, mapq: Int = 60): SamRecord = {
    require(r == 0 || r == 1 || r == 2)
    val rec = builder.addFrag(contig=builder.dict(chrom).index, start=pos, strand=strand, cigar=cigar, mapq=mapq).get
    rec.supplementary = supp

    if (r > 0) {
      rec.paired = true
      rec.firstOfPair  = r == 1
      rec.secondOfPair = r == 2
    }

    rec
  }

  /** Construct a Template from one or more SamRecords. */
  def t(recs: SamRecord*): Template = Template(recs.iterator.tapEach(_.name = "q1"))

  /** Call breakpoints from a template. */
  def call(t: Template): IndexedSeq[BreakpointEvidence] =
    SvPileup.findBreakpoints(
      template                       = t,
      maxWithinReadDistance          = 5,
      maxReadPairInnerDistance       = 1000,
      minUniqueBasesToAdd            = 10
    )

  /** Short hand for constructing a BreakpointEvidence. */
  def bp(ev: EvidenceType,
         lChrom: String, lPos: Int, lStrand: Strand,
         rChrom: String, rPos: Int, rStrand: Strand,
         recs: Iterator[SamRecord] = Iterator.empty): BreakpointEvidence = {
    BreakpointEvidence(
      Breakpoint(
        leftRefIndex  = builder.dict(lChrom).index,
        leftPos       = lPos,
        leftPositive  = lStrand == Plus,
        rightRefIndex = builder.dict(rChrom).index,
        rightPos      = rPos,
        rightPositive = rStrand == Plus),
      ev,
      recs.toSet)
  }


  "SvPileup.findBreakpoints(template)" should "find nothing interesting in an FR read pair with no supplementaries" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="100M", supp=false),
      r("chr1", 250, Minus, r=2, cigar="100M", supp=false),
    )
    call(template) should contain theSameElementsInOrderAs IndexedSeq.empty
  }

  it should "not call a breakpoint from an FR pair with overlapping reads" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="100M", supp=false),
      r("chr1", 150, Minus, r=2, cigar="100M", supp=false),
    )
    call(template) should contain theSameElementsInOrderAs IndexedSeq.empty
  }

  it should "call a breakpoint from a tandem read pair" in {
    val template = t(
      r("chr1", 100, Plus, r=1, cigar="100M", supp=false),
      r("chr1", 250, Plus, r=2, cigar="100M", supp=false),
    )
    val bps = call(template)

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 199, Plus, "chr1", 349, Minus, template.allReads)
    )
  }

  it should "call a breakpoint from an RF read pair" in {
    val template = t(
      r("chr1", 100, Minus, r=1, cigar="100M", supp=false),
      r("chr1", 250, Plus,  r=2, cigar="100M", supp=false),
    )
    val bps = call(template)

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 100, Minus, "chr1", 349, Minus, template.allReads)
    )
  }

  it should "call a breakpoint from an FR read pair with a large insert size" in {
    val template = t(
      r("chr1", 100,   Plus,   r=1, cigar="100M", supp=false),
      r("chr1", 10000, Minus,  r=2, cigar="100M", supp=false),
    )
    val bps = call(template)

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 199, Plus, "chr1", 10000, Plus, template.allReads)
    )
  }

  it should "call a breakpoint from an FR read pair across chromosomes" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="100M", supp=false),
      r("chr2", 300, Minus, r=2, cigar="100M", supp=false),
    )
    val bps = call(template)

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 199, Plus, "chr2", 300, Plus, template.allReads)
    )
  }

  it should "call a breakpoint with one or more split reads" in {
    // A set of reads where each read (r1/r2) is either completely aligned near one side of
    // breakpoint, or the read is split-read aligned across the breakpoint.
    val fullR1  = r("chr1",   1, Plus,  r=1, cigar="100M",   supp=false)
    val r1Half1 = r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false)
    val r1Half2 = r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true)

    val fullR2  = r("chr7", 850, Minus, r=2, cigar="100M",   supp=false)
    val r2Half1 = r("chr7", 800, Minus, r=2, cigar="30S70M", supp=false)
    val r2Half2 = r("chr1", 120, Minus, r=2, cigar="30M70S", supp=true)

    call(t(r1Half2, r1Half1, fullR2)) should contain theSameElementsInOrderAs IndexedSeq(
      bp(SplitRead, "chr1", 149, Plus, "chr7", 800, Plus, Iterator(r1Half1, r1Half2))
    )
    call(t(fullR1,  r2Half1, r2Half2)) should contain theSameElementsInOrderAs IndexedSeq(
      bp(SplitRead, "chr1", 149, Plus, "chr7", 800, Plus, Iterator(r2Half1, r2Half2))
    )
    call(t(r1Half1, r1Half2, r2Half1, r2Half2)) should contain theSameElementsInOrderAs IndexedSeq(
      bp(SplitRead, "chr1", 149, Plus, "chr7", 800, Plus, Iterator(r1Half1, r1Half2, r2Half1, r2Half2))
    )
  }

  it should "only call one breakpoint where R1 and R2 have slightly different split-read support" in {
    // If there is repetitive sequence at the breakpoint the aligner split r1 and r2 slightly
    // differently on each side - make sure we don't generate two breakpoints in that case
    val template = t(
      // Read 1 supports a breakpoint at chr1:149F>chr7:800F
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true),
      // Read 1 supports a breakpoint at chr1:151F>chr7:802F
      r("chr7", 802, Minus, r=2, cigar="30S70M", supp=false),
      r("chr1", 122, Minus, r=2, cigar="30M70S", supp=true)
    )
    val bps = call(template)

    // The one breakpoint that does come out is not ideal in this case, but maybe that's ok?
    // I think ideally it would call chr1:149F>chr7:800F or chr1:151F>chr7:802F not chr1:151F>chr7:800F
    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(SplitRead, "chr1", 151, Plus, "chr7", 800, Plus, template.allReads)
    )
  }

  it should "call more than one breakpoint from a read pair that supports more than one" in {
    val recs = IndexedSeq(
      r("chr1", 100, Plus,  r=1, cigar="30M70S",    supp=false),
      r("chr2", 500, Minus, r=1, cigar="30S40M30S", supp=true),
      r("chr3", 900, Plus,  r=1, cigar="70S30M",    supp=true),

      r("chr3", 1200, Minus, r=2, cigar="100M", supp=false),
    )
    val bps = call(t(recs:_*))

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(SplitRead, "chr1", 129, Plus,  "chr2", 539, Minus, recs.take(2).iterator),
      bp(SplitRead, "chr2", 500, Minus, "chr3", 900, Plus, recs.slice(1, 3).iterator),
    )
  }

  it should "call a breakpoint from a single-end split read with no mate" in {
    val r1Half1 = r("chr1", 100, Plus,  r=0, cigar="50M50S", supp=false)
    val r1Half2 = r("chr7", 800, Plus,  r=0, cigar="50S50M", supp=true)

    val expected = IndexedSeq(bp(SplitRead, "chr1", 149, Plus, "chr7", 800, Plus, Iterator(r1Half1, r1Half2)))
    call(t(r1Half1, r1Half2)) should contain theSameElementsInOrderAs expected
  }

  "SvPileup.filterTemplate" should "do nothing to a template with just high mapq primaries" in {
    val template = t(
      r("chr1", 100, Plus, r=1, cigar="100M", supp=false, mapq=50),
      r("chr1", 300, Plus, r=2, cigar="100M", supp=false, mapq=50),
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe template
  }

  it should "do nothing to a template with high mapq primaries and supplementaries" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false, mapq=50),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true,  mapq=50),
      r("chr7", 800, Minus, r=2, cigar="30S70M", supp=false, mapq=50),
      r("chr1", 120, Minus, r=2, cigar="30M70S", supp=true,  mapq=50),
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe template
  }

  it should "remove a low quality supplementary record" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false, mapq=50),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true,  mapq=50),
      r("chr7", 800, Minus, r=2, cigar="30S70M", supp=false, mapq=50),
      r("chr1", 120, Minus, r=2, cigar="30M70S", supp=true,  mapq=1),
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe template.copy(r2Supplementals=Nil)
  }

  it should "remove all evidence of R2 if the primary mapping is low quality" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false, mapq=50),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true,  mapq=50),
      r("chr7", 800, Minus, r=2, cigar="30S70M", supp=false, mapq=1),
      r("chr1", 120, Minus, r=2, cigar="30M70S", supp=true,  mapq=50),
    )

    val expected = template.copy(r2=None, r2Supplementals=Nil)
    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe expected
  }

  it should "return None if both R1 and R2 primaries are low quality" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false, mapq=1),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true,  mapq=50),
      r("chr7", 800, Minus, r=2, cigar="30S70M", supp=false, mapq=1),
      r("chr1", 120, Minus, r=2, cigar="30M70S", supp=true,  mapq=50),
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20) shouldBe None
  }

  it should "remove an unmapped R2" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false, mapq=50),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true,  mapq=50),
      Seq(r("chr1", 100, Plus, r=2, cigar="100M",   supp=false, mapq=0)).tapEach(_.unmapped = true).head,
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe template.copy(r2=None)
  }

  it should "handle a template with fragment data with high mapping quality" in {
    val template = t(
      r("chr1", 100, Plus,  r=0, cigar="50M50S", supp=false, mapq=50),
      r("chr7", 800, Plus,  r=0, cigar="50S50M", supp=true,  mapq=50),
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe template
  }

  it should "remove low mapq supplementary records from a fragment template" in {
    val template = t(
      r("chr1", 100, Plus,  r=0, cigar="50M50S",    supp=false, mapq=50),
      r("chr7", 800, Minus, r=0, cigar="50S20M30S", supp=true,  mapq=7),
      r("chr7", 820, Plus,  r=0, cigar="70S30M",    supp=true,  mapq=50),
    )

    val expected = template.copy(r1Supplementals = template.r1Supplementals.filter(_.start != 800))
    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20).value shouldBe expected
  }

  it should "return None for a fragment template with a low quality primary mapping" in {
    val template = t(
      r("chr1", 100, Plus,  r=0, cigar="50M50S",    supp=false, mapq=5),
      r("chr7", 800, Minus, r=0, cigar="50S20M30S", supp=true,  mapq=50),
      r("chr7", 820, Plus,  r=0, cigar="70S30M",    supp=true,  mapq=50),
    )

    SvPileup.filterTemplate(template, minPrimaryMapq=30, minSupplementaryMapq=20) shouldBe None
  }

  /** Writes the given BAM records to a BAM file */
  private def toInput(rec: SamRecord*): PathToBam = {
    val path    = makeTempFile("input.", ".bam")
    val header  = new SamBuilder().header
    val writer  = SamWriter(path, header, sort=None)
    writer ++= rec
    writer.close()
    path
  }

  private object Outputs {
    def apply(): Outputs = {
      val prefix = makeTempFile("output", "")
      Outputs(
        prefix = prefix,
        bam    = PathUtil.pathTo(prefix + ".bam"),
        txt    = PathUtil.pathTo(prefix + ".bam")
      )
    }
  }

  private case class Outputs(prefix: PathPrefix, bam: PathToBam, txt: FilePath) {
    Seq(prefix, bam, txt)foreach(_.toFile.deleteOnExit())
  }

  "SvPileup" should "run end to end" in {
    // A set of reads where each read (r1/r2) is either completely aligned near one side of
    // breakpoint, or the read is split-read aligned across the breakpoint.
    val fullR1  = r("chr1",   1, Plus,  r=1, cigar="100M",   supp=false)
    val r1Half1 = r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false)
    val r1Half2 = r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true)

    val fullR2  = r("chr7", 850, Minus, r=2, cigar="100M",   supp=false)
    val r2Half1 = r("chr7", 800, Minus, r=2, cigar="30S70M", supp=false)
    val r2Half2 = r("chr1", 120, Minus, r=2, cigar="30M70S", supp=true)

    Seq(fullR1, r1Half1, r1Half2, fullR2, r2Half1, r2Half2).foreach { rec => rec.name = "q1" }

    def test(recs: Seq[SamRecord], ids: Seq[String]): Unit = {
      recs.length shouldBe ids.length
      val input   = toInput(recs:_*)
      val outputs = Outputs()
      new SvPileup(input=input, output=outputs.prefix).execute()
      val outRecs = readBamRecs(outputs.bam)
      outRecs.length shouldBe recs.length
      outRecs.zip(ids).foreach { case (rec, id) =>
        withClue(rec.asSam.getSAMString.trim) {
          rec.getOrElse[String](SvPileup.SamBreakpointTag, "") shouldBe id
        }
      }
    }

    // r1Half2 is not annotated, since it is superseded by r1Half1 -> fullR2
    SamPairUtil.setMateInfo(r1Half1.asSam, fullR2.asSam, true)
    SamPairUtil.setMateInformationOnSupplementalAlignment(r1Half2.asSam, fullR2.asSam, true)
    test(Seq(r1Half2, r1Half1, fullR2), Seq("0", "", "0"))

    // fullR1 does not contain the breakpoint, while r2Half1 -> r2Half2 does
    SamPairUtil.setMateInfo(fullR1.asSam, r2Half1.asSam, true)
    SamPairUtil.setMateInformationOnSupplementalAlignment(r2Half2.asSam, fullR1.asSam, true)
    test(Seq(fullR1,  r2Half1, r2Half2), Seq("", "0", "0"))

    // all the reads are annotated!
    SamPairUtil.setMateInfo(r1Half1.asSam, r2Half1.asSam, true)
    SamPairUtil.setMateInformationOnSupplementalAlignment(r2Half2.asSam, r1Half1.asSam, true)
    SamPairUtil.setMateInformationOnSupplementalAlignment(r1Half2.asSam, r2Half1.asSam, true)
    test(Seq(r1Half1, r1Half2, r2Half1, r2Half2), Seq("0", "0", "0", "0"))
  }
}
