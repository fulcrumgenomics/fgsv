package com.fulcrumgenomics.fgsv.tools

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.fgsv.UnitSpec
import com.fulcrumgenomics.fgsv.util.BlockOrigin.{Both, ReadOne, ReadTwo}
import com.fulcrumgenomics.fgsv.util.EvidenceType._
import com.fulcrumgenomics.fgsv.util.{AlignmentBlock, BlockOrigin, GenomicRange, PutativeBreakpoint}

class SVPileupTest extends UnitSpec {
  private def fromRangeOnly(refIndex: Int, start: Int, end: Int, origin: BlockOrigin = ReadOne): AlignmentBlock = AlignmentBlock(
    origin=origin,
    start=1, end=100, positiveStrand=true, cigar=Cigar("100M"),
    range=GenomicRange(refIndex=refIndex, start=start, end=end)
  )

  "SVPileup.determineInterContigChimera" should "determine if two alignment blocks are an inter-contig chimera" in {
    val chr1 = fromRangeOnly(refIndex=0, start=1, end=100)
    val chr2 = fromRangeOnly(refIndex=1, start=1, end=100)

    // Same contig
    SVPileup.determineInterContigChimera(b1=chr1, b2=chr1).isEmpty shouldBe true

    // Different contig
    SVPileup.determineInterContigChimera(b1=chr1, b2=chr2).value shouldBe PutativeBreakpoint(
      chr1, chr2, SplitReadInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr2, b2=chr1).value shouldBe PutativeBreakpoint(
      chr2, chr1, SplitReadInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr1, b2=chr2.copy(origin=ReadTwo)).value shouldBe PutativeBreakpoint(
      chr1, chr2.copy(origin=ReadTwo), ReadPairInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr2.copy(origin=ReadTwo), b2=chr1).value shouldBe PutativeBreakpoint(
      chr2.copy(origin=ReadTwo), chr1, ReadPairInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr1, b2=chr2.copy(origin=Both)).value shouldBe PutativeBreakpoint(
      chr1, chr2.copy(origin=Both), ReadPairInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr2.copy(origin=Both), b2=chr1).value shouldBe PutativeBreakpoint(
      chr2.copy(origin=Both), chr1, ReadPairInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr1.copy(origin=Both), b2=chr2.copy(origin=Both)).value shouldBe PutativeBreakpoint(
      chr1.copy(origin=Both), chr2.copy(origin=Both), ReadPairInterContig
    )
    SVPileup.determineInterContigChimera(b1=chr2.copy(origin=Both), b2=chr1.copy(origin=Both)).value shouldBe PutativeBreakpoint(
      chr2.copy(origin=Both), chr1.copy(origin=Both), ReadPairInterContig
    )
  }

  "SVPileup.determineOddPairOrientation" should "determine if two alignment blocks have an odd pair orientation" in {
    val earlier = fromRangeOnly(refIndex=1, start=1, end=100)
    val overlap = fromRangeOnly(refIndex=1, start=50, end=150)
    val later   = fromRangeOnly(refIndex=1, start=150, end=200)

    // same block
    SVPileup.determineIntraContigChimera(b1=earlier, b2=earlier, maxInnerDistance=0).isEmpty shouldBe true

    // overlapping blocks
    SVPileup.determineIntraContigChimera(b1=earlier, b2=overlap, maxInnerDistance=0).isEmpty shouldBe true

    // non-overlapping blocks
    SVPileup.determineIntraContigChimera(b1=earlier, b2=later, maxInnerDistance=49).value shouldBe PutativeBreakpoint(
      earlier, later, SplitReadIntraContig
    )
    SVPileup.determineIntraContigChimera(b1=earlier, b2=later.copy(origin=ReadTwo), maxInnerDistance=49).value shouldBe PutativeBreakpoint(
      earlier, later.copy(origin=ReadTwo), ReadPairIntraContig
    )
    SVPileup.determineIntraContigChimera(b1=earlier, b2=later.copy(origin=Both), maxInnerDistance=49).value shouldBe PutativeBreakpoint(
      earlier, later.copy(origin=Both), ReadPairIntraContig
    )
    SVPileup.determineIntraContigChimera(b1=earlier.copy(origin=Both), b2=later, maxInnerDistance=49).value shouldBe PutativeBreakpoint(
      earlier.copy(origin=Both), later, ReadPairIntraContig
    )
    SVPileup.determineIntraContigChimera(b1=earlier.copy(origin=Both), b2=later.copy(origin=Both), maxInnerDistance=49).value shouldBe PutativeBreakpoint(
      earlier.copy(origin=Both), later.copy(origin=Both), ReadPairIntraContig
    )
    SVPileup.determineIntraContigChimera(b1=earlier, b2=later, maxInnerDistance=50).isEmpty shouldBe true
  }

  "SVPileup.determineOddPairOrientation" should "determine if two alignment blocks have an odd read-pair orientation" in {
    val r1Forward = fromRangeOnly(refIndex=1, start=1, end=100).copy(positiveStrand=true)
    val r1Reverse = fromRangeOnly(refIndex=1, start=100, end=200).copy(positiveStrand=false)
    val r2Forward = fromRangeOnly(refIndex=1, start=1, end=100).copy(positiveStrand=true, origin=ReadTwo)
    val r2Reverse = fromRangeOnly(refIndex=1, start=100, end=200).copy(positiveStrand=false, origin=ReadTwo)
    val bothForward = fromRangeOnly(refIndex=1, start=1, end=100).copy(positiveStrand=true, origin=Both)
    val bothReverse = fromRangeOnly(refIndex=1, start=1, end=200).copy(positiveStrand=false, origin=Both)


    // read pairs expect to be FR, so no breakpoint
    SVPileup.determineOddPairOrientation(b1=r1Forward, b2=r2Reverse).isEmpty shouldBe true
    SVPileup.determineOddPairOrientation(b1=r2Forward, b2=r1Reverse).isEmpty shouldBe true
    SVPileup.determineOddPairOrientation(b1=r1Reverse, b2=r2Forward).isEmpty shouldBe true
    SVPileup.determineOddPairOrientation(b1=r2Reverse, b2=r1Forward).isEmpty shouldBe true

    // split pairs expect to be FF or RR, so no breakpoint
    SVPileup.determineOddPairOrientation(b1=r1Forward, b2=r1Forward).isEmpty shouldBe true
    SVPileup.determineOddPairOrientation(b1=r1Reverse, b2=r1Reverse).isEmpty shouldBe true

    // read pairs should be odd if RF or tandem
    SVPileup.determineOddPairOrientation(b1=r1Forward, b2=r2Forward).value shouldBe PutativeBreakpoint(
      r1Forward, r2Forward, ReadPairTandem
    )
    SVPileup.determineOddPairOrientation(b1=r1Reverse, b2=r2Reverse).value shouldBe PutativeBreakpoint(
      r1Reverse, r2Reverse, ReadPairTandem
    )
    // NB: bothForward can be treated as r2Forward, so tandem
    SVPileup.determineOddPairOrientation(b1=r1Forward, b2=bothForward).value shouldBe PutativeBreakpoint(
      r1Forward, bothForward, ReadPairTandem
    )
    // NB: bothReverse can be treated as r2Reverse, so tandem
    SVPileup.determineOddPairOrientation(b1=r1Reverse, b2=bothReverse).value shouldBe PutativeBreakpoint(
      r1Reverse, bothReverse, ReadPairTandem
    )
    // NB: both* can be treated as R1/R2, so tandem
    SVPileup.determineOddPairOrientation(b1=bothForward, b2=bothForward).value shouldBe PutativeBreakpoint(
      bothForward, bothForward, ReadPairTandem
    )
    // NB: both* can be treated as R1/R2, so tandem
    SVPileup.determineOddPairOrientation(b1=bothReverse, b2=bothReverse).value shouldBe PutativeBreakpoint(
      bothReverse, bothReverse, ReadPairTandem
    )

    // split reads should be odd if RF or FR
    Seq(
      (r1Forward, r1Reverse),
      (r1Reverse, r1Forward),
      (r2Forward, r2Reverse),
      (r2Reverse, r2Forward),
    ).foreach { case (b1, b2) =>
      SVPileup.determineOddPairOrientation(b1=b1, b2=b2).value shouldBe PutativeBreakpoint(
        b1, b2, SplitReadOppositeStrand
      )
    }
  }
}
