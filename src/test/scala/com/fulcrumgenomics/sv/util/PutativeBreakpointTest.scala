package com.fulcrumgenomics.sv.util

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.sv.UnitSpec

class PutativeBreakpointTest extends UnitSpec {
  // Constructs an alignment block that is all-matching with the given starts and length
  private def block(readStart: Int, refIndex: Int, refStart: Int, length: Int, positive: Boolean): AlignedSegment = {
    AlignedSegment(
      origin         = SegmentOrigin.ReadOne,
      readStart      = readStart,
      readEnd        = readStart + length - 1,
      positiveStrand = positive,
      cigar          = Cigar(s"${length}M"),
      range          = GenomicRange(refIndex, refStart, refStart + length - 1)
    )
  }

  "PutativeBreakpoint.apply(block, block)" should "correctly capture a simple F-F breakpoint" in {
    // F-F breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val block1 = block(readStart=1,  refIndex=1, refStart=100,  length=50, positive=true)
    val block2 = block(readStart=51, refIndex=1, refStart=1000, length=30, positive=true)
    val bp     = PutativeBreakpoint(block1, block2, EvidenceType.SplitReadIntraContig)
    (bp.leftPos, bp.leftPositive)   shouldBe (149, true)
    (bp.rightPos, bp.rightPositive) shouldBe (1000, true)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revBlock1 = block(readStart=1,  refIndex=1, refStart=1000, length=30, positive=false)
    val revBlock2 = block(readStart=31, refIndex=1, refStart=100,  length=50, positive=false)
    val rbp       = PutativeBreakpoint(revBlock1, revBlock2, EvidenceType.SplitReadIntraContig)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (149, true)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1000, true)
  }

  it should "correctly capture an F-R breakpoint" in {
    // F-R breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val block1 = block(readStart=1,  refIndex=1, refStart=100,  length=50, positive=true)
    val block2 = block(readStart=51, refIndex=1, refStart=1000, length=30, positive=false)
    val bp     = PutativeBreakpoint(block1, block2, EvidenceType.SplitReadIntraContig)
    (bp.leftPos, bp.leftPositive)   shouldBe (149, true)
    (bp.rightPos, bp.rightPositive) shouldBe (1029, false)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revBlock1 = block(readStart=1, refIndex=1, refStart=1000, length=30, positive=true)
    val revBlock2 = block(readStart=31,  refIndex=1, refStart=100,  length=50, positive=false)
    val rbp       = PutativeBreakpoint(revBlock1, revBlock2, EvidenceType.SplitReadIntraContig)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (149, true)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1029, false)
  }

  it should "correctly capture an R-R breakpoint" in {
    // R-R breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val block1 = block(readStart=1,  refIndex=1, refStart=100,  length=50, positive=false)
    val block2 = block(readStart=51, refIndex=1, refStart=1000, length=30, positive=false)
    val bp     = PutativeBreakpoint(block1, block2, EvidenceType.SplitReadIntraContig)
    (bp.leftPos, bp.leftPositive)   shouldBe (100, false)
    (bp.rightPos, bp.rightPositive) shouldBe (1029, false)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revBlock1 = block(readStart=1, refIndex=1, refStart=1000, length=30, positive=true)
    val revBlock2 = block(readStart=31,  refIndex=1, refStart=100,  length=50, positive=true)
    val rbp       = PutativeBreakpoint(revBlock1, revBlock2, EvidenceType.SplitReadIntraContig)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (100, false)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1029, false)

  }

  it should "correctly capture an R-F breakpoint" in {
    // R-F breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val block1 = block(readStart=1,  refIndex=1, refStart=100,  length=50, positive=false)
    val block2 = block(readStart=51, refIndex=1, refStart=1000, length=30, positive=true)
    val bp     = PutativeBreakpoint(block1, block2, EvidenceType.SplitReadIntraContig)
    (bp.leftPos, bp.leftPositive)   shouldBe (100, false)
    (bp.rightPos, bp.rightPositive) shouldBe (1000, true)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revBlock1 = block(readStart=1, refIndex=1, refStart=1000, length=30, positive=false)
    val revBlock2 = block(readStart=31,  refIndex=1, refStart=100,  length=50, positive=true)
    val rbp       = PutativeBreakpoint(revBlock1, revBlock2, EvidenceType.SplitReadIntraContig)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (100, false)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1000, true)
  }

  it should "represent a tandem duplication sanely" in {
    // We want the same representation regardless of which piece of the duplication was sequenced more so
    // we test where the second block is shorter, the same length as, and longer than the first
    Seq(-1, 0, 1).foreach { addend =>
      val block1Length = 50
      val block2Length = block1Length + addend
      val block1 = block(readStart=1,  refIndex=1, refStart=100, length=block1Length, positive=true)
      val block2 = block(readStart=51, refIndex=1, refStart=100, length=block2Length, positive=true)
      val bp     = PutativeBreakpoint(block1, block2, EvidenceType.SplitReadIntraContig)
      (bp.leftPos,  bp.leftPositive)  shouldBe (100, false)
      (bp.rightPos, bp.rightPositive) shouldBe (149, false)
    }
  }
}