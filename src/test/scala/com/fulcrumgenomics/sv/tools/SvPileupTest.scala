package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.sv.EvidenceType._
import com.fulcrumgenomics.sv.SegmentOrigin.{Both, ReadOne, ReadTwo}
import com.fulcrumgenomics.sv._

class SvPileupTest extends UnitSpec {
  private def fromRangeOnly(refIndex: Int, start: Int, end: Int, origin: SegmentOrigin = ReadOne): AlignedSegment = AlignedSegment(
    origin=origin,
    readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar("100M"),
    range=GenomicRange(refIndex=refIndex, start=start, end=end)
  )

  "SvPileup.determineInterContigBreakpiont" should "determine if two alignment segments are an inter-contig breakpoint" in {
    val chr1 = fromRangeOnly(refIndex=0, start=1, end=100)
    val chr2 = fromRangeOnly(refIndex=1, start=1, end=100)

    // Same contig
    SvPileup.determineInterContigBreakpiont(seg1=chr1, seg2=chr1).isEmpty shouldBe true

    // Different contig
    SvPileup.determineInterContigBreakpiont(seg1=chr1, seg2=chr2).value shouldBe BreakpointEvidence(
      chr1, chr2, SplitReadInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr2, seg2=chr1).value shouldBe BreakpointEvidence(
      chr2, chr1, SplitReadInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr1, seg2=chr2.copy(origin=ReadTwo)).value shouldBe BreakpointEvidence(
      chr1, chr2.copy(origin=ReadTwo), ReadPairInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr2.copy(origin=ReadTwo), seg2=chr1).value shouldBe BreakpointEvidence(
      chr2.copy(origin=ReadTwo), chr1, ReadPairInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr1, seg2=chr2.copy(origin=Both)).value shouldBe BreakpointEvidence(
      chr1, chr2.copy(origin=Both), ReadPairInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr2.copy(origin=Both), seg2=chr1).value shouldBe BreakpointEvidence(
      chr2.copy(origin=Both), chr1, ReadPairInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr1.copy(origin=Both), seg2=chr2.copy(origin=Both)).value shouldBe BreakpointEvidence(
      chr1.copy(origin=Both), chr2.copy(origin=Both), ReadPairInterContig
    )
    SvPileup.determineInterContigBreakpiont(seg1=chr2.copy(origin=Both), seg2=chr1.copy(origin=Both)).value shouldBe BreakpointEvidence(
      chr2.copy(origin=Both), chr1.copy(origin=Both), ReadPairInterContig
    )
  }

  "SvPileup.determineOddPairOrientation" should "determine if two alignment segments have an odd pair orientation" in {
    val earlier = fromRangeOnly(refIndex=1, start=1, end=100)
    val overlap = fromRangeOnly(refIndex=1, start=50, end=150)
    val later   = fromRangeOnly(refIndex=1, start=150, end=200)

    // same segment
    SvPileup.determineIntraContigBreakpiont(seg1=earlier, seg2=earlier, maxInnerDistance=0).isEmpty shouldBe true

    // overlapping segments
    SvPileup.determineIntraContigBreakpiont(seg1=earlier, seg2=overlap, maxInnerDistance=0).isEmpty shouldBe true

    // non-overlapping segments
    SvPileup.determineIntraContigBreakpiont(seg1=earlier, seg2=later, maxInnerDistance=49).value shouldBe BreakpointEvidence(
      earlier, later, SplitReadIntraContig
    )
    SvPileup.determineIntraContigBreakpiont(seg1=earlier, seg2=later.copy(origin=ReadTwo), maxInnerDistance=49).value shouldBe BreakpointEvidence(
      earlier, later.copy(origin=ReadTwo), ReadPairIntraContig
    )
    SvPileup.determineIntraContigBreakpiont(seg1=earlier, seg2=later.copy(origin=Both), maxInnerDistance=49).value shouldBe BreakpointEvidence(
      earlier, later.copy(origin=Both), ReadPairIntraContig
    )
    SvPileup.determineIntraContigBreakpiont(seg1=earlier.copy(origin=Both), seg2=later, maxInnerDistance=49).value shouldBe BreakpointEvidence(
      earlier.copy(origin=Both), later, ReadPairIntraContig
    )
    SvPileup.determineIntraContigBreakpiont(seg1=earlier.copy(origin=Both), seg2=later.copy(origin=Both), maxInnerDistance=49).value shouldBe BreakpointEvidence(
      earlier.copy(origin=Both), later.copy(origin=Both), ReadPairIntraContig
    )
    SvPileup.determineIntraContigBreakpiont(seg1=earlier, seg2=later, maxInnerDistance=50).isEmpty shouldBe true
  }

  "SvPileup.determineOddPairOrientation" should "determine if two alignment segments have an odd read-pair orientation" in {
    val r1Forward = fromRangeOnly(refIndex=1, start=1, end=100).copy(positiveStrand=true)
    val r1Reverse = fromRangeOnly(refIndex=1, start=100, end=200).copy(positiveStrand=false)
    val r2Forward = fromRangeOnly(refIndex=1, start=1, end=100).copy(positiveStrand=true, origin=ReadTwo)
    val r2Reverse = fromRangeOnly(refIndex=1, start=100, end=200).copy(positiveStrand=false, origin=ReadTwo)
    val bothForward = fromRangeOnly(refIndex=1, start=1, end=100).copy(positiveStrand=true, origin=Both)
    val bothReverse = fromRangeOnly(refIndex=1, start=1, end=200).copy(positiveStrand=false, origin=Both)


    // read pairs expect to be FR, so no breakpoint
    SvPileup.determineOddPairOrientation(seg1=r1Forward, seg2=r2Reverse).isEmpty shouldBe true
    SvPileup.determineOddPairOrientation(seg1=r2Forward, seg2=r1Reverse).isEmpty shouldBe true
    SvPileup.determineOddPairOrientation(seg1=r1Reverse, seg2=r2Forward).isEmpty shouldBe true
    SvPileup.determineOddPairOrientation(seg1=r2Reverse, seg2=r1Forward).isEmpty shouldBe true

    // split pairs expect to be FF or RR, so no breakpoint
    SvPileup.determineOddPairOrientation(seg1=r1Forward, seg2=r1Forward).isEmpty shouldBe true
    SvPileup.determineOddPairOrientation(seg1=r1Reverse, seg2=r1Reverse).isEmpty shouldBe true

    // read pairs should be odd if RF or tandem
    SvPileup.determineOddPairOrientation(seg1=r1Forward, seg2=r2Forward).value shouldBe BreakpointEvidence(
      r1Forward, r2Forward, ReadPairTandem
    )
    SvPileup.determineOddPairOrientation(seg1=r1Reverse, seg2=r2Reverse).value shouldBe BreakpointEvidence(
      r1Reverse, r2Reverse, ReadPairTandem
    )
    // NB: bothForward can be treated as r2Forward, so tandem
    SvPileup.determineOddPairOrientation(seg1=r1Forward, seg2=bothForward).value shouldBe BreakpointEvidence(
      r1Forward, bothForward, ReadPairTandem
    )
    // NB: bothReverse can be treated as r2Reverse, so tandem
    SvPileup.determineOddPairOrientation(seg1=r1Reverse, seg2=bothReverse).value shouldBe BreakpointEvidence(
      r1Reverse, bothReverse, ReadPairTandem
    )
    // NB: both* can be treated as R1/R2, so tandem
    SvPileup.determineOddPairOrientation(seg1=bothForward, seg2=bothForward).value shouldBe BreakpointEvidence(
      bothForward, bothForward, ReadPairTandem
    )
    // NB: both* can be treated as R1/R2, so tandem
    SvPileup.determineOddPairOrientation(seg1=bothReverse, seg2=bothReverse).value shouldBe BreakpointEvidence(
      bothReverse, bothReverse, ReadPairTandem
    )

    // split reads should be odd if RF or FR
    Seq(
      (r1Forward, r1Reverse),
      (r1Reverse, r1Forward),
      (r2Forward, r2Reverse),
      (r2Reverse, r2Forward),
    ).foreach { case (seg1, seg2) =>
      SvPileup.determineOddPairOrientation(seg1=seg1, seg2=seg2).value shouldBe BreakpointEvidence(
        seg1, seg2, SplitReadOppositeStrand
      )
    }
  }
}
