package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.Template
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.sv.EvidenceType.{ReadPair, SplitRead}
import com.fulcrumgenomics.sv.SegmentOrigin.{Both, ReadOne, ReadTwo}
import com.fulcrumgenomics.sv._
import com.fulcrumgenomics.testing.SamBuilder

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
  def r(chrom: String, pos: Int, strand: SamBuilder.Strand, r: Int, cigar: String, supp: Boolean): SamRecord = {
    require(r == 1 || r == 2)
    val rec = builder.addFrag(contig=builder.dict(chrom).index, start=pos, strand=strand, cigar=cigar).get
    rec.supplementary = supp
    rec.paired = true
    rec.firstOfPair = r == 1
    rec.secondOfPair = r == 2
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
      minSupplementaryMappingQuality = 0,
      minUniqueBasesToAdd            = 10
    )

  /** Short hand for constructing a BreakpointEvidence. */
  def bp(ev: EvidenceType, lChrom: String, lPos: Int, lStrand: Strand, rChrom: String, rPos: Int, rStrand: Strand): BreakpointEvidence =
    BreakpointEvidence(
      Breakpoint(
        leftRefIndex  = builder.dict(lChrom).index,
        leftPos       = lPos,
        leftPositive  = lStrand == Plus,
        rightRefIndex = builder.dict(rChrom).index,
        rightPos      = rPos,
        rightPositive = rStrand == Plus),
      ev)



  "SvPileup.findBreakpoints(template)" should "find nothing interesting in an FR read pair with no supplementaries" in {
    val template = t(
      r("chr1", 100, Plus,  r=1, cigar="100M", supp=false),
      r("chr1", 250, Minus, r=2, cigar="100M", supp=false),
    )
    call(template) should contain theSameElementsInOrderAs IndexedSeq.empty
  }

  it should "call a breakpoint from a tandem read pair" in {
    val bps = call(t(
      r("chr1", 100, Plus, r=1, cigar="100M", supp=false),
      r("chr1", 250, Plus, r=2, cigar="100M", supp=false),
    ))

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 199, Plus, "chr1", 349, Minus)
    )
  }

  it should "call a breakpoint from an RF read pair" in {
    val bps = call(t(
      r("chr1", 100, Minus, r=1, cigar="100M", supp=false),
      r("chr1", 250, Plus,  r=2, cigar="100M", supp=false),
    ))

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 100, Minus, "chr1", 349, Minus)
    )
  }

  it should "call a breakpoint from an FR read pair with a large insert size" in {
    val bps = call(t(
      r("chr1", 100,   Plus,   r=1, cigar="100M", supp=false),
      r("chr1", 10000, Minus,  r=2, cigar="100M", supp=false),
    ))

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 199, Plus, "chr1", 10000, Plus)
    )
  }

  it should "call a breakpoint from an FR read pair across chromosomes" in {
    val bps = call(t(
      r("chr1", 100, Plus,  r=1, cigar="100M", supp=false),
      r("chr2", 300, Minus, r=2, cigar="100M", supp=false),
    ))

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(ReadPair, "chr1", 199, Plus, "chr2", 300, Plus)
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

    val expected = IndexedSeq(bp(SplitRead, "chr1", 149, Plus, "chr7", 800, Plus))
    call(t(r1Half2, r1Half1, fullR2))           should contain theSameElementsInOrderAs expected
    call(t(fullR1,  r2Half1, r2Half2))          should contain theSameElementsInOrderAs expected
    call(t(r1Half1, r1Half2, r2Half1, r2Half2)) should contain theSameElementsInOrderAs expected
  }

  it should "only call one breakpoint where R1 and R2 have slightly different split-read support" in {
    // If there is repetitive sequence at the breakpoint the aligner split r1 and r2 slightly
    // differently on each side - make sure we don't generate two breakpoints in that case
    val bps = call(t(
      // Read 1 supports a breakpoint at chr1:149F>chr7:800F
      r("chr1", 100, Plus,  r=1, cigar="50M50S", supp=false),
      r("chr7", 800, Plus,  r=1, cigar="50S50M", supp=true),
      // Read 1 supports a breakpoint at chr1:151F>chr7:802F
      r("chr7", 802, Minus, r=2, cigar="30S70M", supp=false),
      r("chr1", 122, Minus, r=2, cigar="30M70S", supp=true)
    ))

    // The one breakpoint that does come out is not ideal in this case, but maybe that's ok?
    // I think ideally it would call chr1:149F>chr7:800F or chr1:151F>chr7:802F not chr1:151F>chr7:800F
    val expected = IndexedSeq(bp(SplitRead, "chr1", 151, Plus, "chr7", 800, Plus))
    bps should contain theSameElementsInOrderAs expected
  }

  it should "call more than one breakpoint from a read pair that supports more than one" in {
    val bps = call(t(
      r("chr1", 100, Plus,  r=1, cigar="30M70S",    supp=false),
      r("chr2", 500, Minus, r=1, cigar="30S40M30S", supp=true),
      r("chr3", 900, Plus,  r=1, cigar="70S30M",    supp=true),

      r("chr3", 1200, Minus, r=2, cigar="100M", supp=false),
    ))

    bps should contain theSameElementsInOrderAs IndexedSeq(
      bp(SplitRead, "chr1", 129, Plus,  "chr2", 539, Minus),
      bp(SplitRead, "chr2", 500, Minus, "chr3", 900, Plus),
    )
  }
}
