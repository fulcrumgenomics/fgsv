package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.sv.SegmentOrigin.{Both, ReadOne, ReadTwo}
import com.fulcrumgenomics.sv._

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
}
