package com.fulcrumgenomics.sv

import com.fulcrumgenomics.alignment.Cigar

class BreakpointTest extends UnitSpec {
  // Constructs an alignment segment that is all-matching with the given starts and length
  private def seg(readStart: Int, refIndex: Int, refStart: Int, length: Int, positive: Boolean): AlignedSegment = {
    AlignedSegment(
      origin         = SegmentOrigin.ReadOne,
      readStart      = readStart,
      readEnd        = readStart + length - 1,
      positiveStrand = positive,
      cigar          = Cigar(s"${length}M"),
      range          = GenomicRange(refIndex, refStart, refStart + length - 1)
    )
  }

  "Breakpoint.apply(segment, segment)" should "correctly capture a simple F-F breakpoint" in {
    // F-F breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val segment1 = seg(readStart=1,  refIndex=1, refStart=100,  length=50, positive=true)
    val segment2 = seg(readStart=51, refIndex=1, refStart=1000, length=30, positive=true)
    val bp       = Breakpoint(segment1, segment2)
    (bp.leftPos, bp.leftPositive)   shouldBe (149, true)
    (bp.rightPos, bp.rightPositive) shouldBe (1000, true)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revSegment1 = seg(readStart=1,  refIndex=1, refStart=1000, length=30, positive=false)
    val revSegment2 = seg(readStart=31, refIndex=1, refStart=100,  length=50, positive=false)
    val rbp         = Breakpoint(revSegment1, revSegment2)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (149, true)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1000, true)
  }

  it should "correctly capture an F-R breakpoint" in {
    // F-R breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val segment1 = seg(readStart=1,  refIndex=1, refStart=100,  length=50, positive=true)
    val segment2 = seg(readStart=51, refIndex=1, refStart=1000, length=30, positive=false)
    val bp       = Breakpoint(segment1, segment2)
    (bp.leftPos, bp.leftPositive)   shouldBe (149, true)
    (bp.rightPos, bp.rightPositive) shouldBe (1029, false)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revSegment1 = seg(readStart=1, refIndex=1, refStart=1000, length=30, positive=true)
    val revSegment2 = seg(readStart=31,  refIndex=1, refStart=100,  length=50, positive=false)
    val rbp         = Breakpoint(revSegment1, revSegment2)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (149, true)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1029, false)
  }

  it should "correctly capture an R-R breakpoint" in {
    // R-R breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val segment1 = seg(readStart=1,  refIndex=1, refStart=100,  length=50, positive=false)
    val segment2 = seg(readStart=51, refIndex=1, refStart=1000, length=30, positive=false)
    val bp       = Breakpoint(segment1, segment2)
    (bp.leftPos, bp.leftPositive)   shouldBe (100, false)
    (bp.rightPos, bp.rightPositive) shouldBe (1029, false)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revSegment1 = seg(readStart=1, refIndex=1, refStart=1000, length=30, positive=true)
    val revSegment2 = seg(readStart=31,  refIndex=1, refStart=100,  length=50, positive=true)
    val rbp         = Breakpoint(revSegment1, revSegment2)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (100, false)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1029, false)
  }

  it should "correctly capture an R-F breakpoint" in {
    // R-F breakpoint if observed with a "forward" read where position in read is correlated to position in genome
    val segment1 = seg(readStart=1,  refIndex=1, refStart=100,  length=50, positive=false)
    val segment2 = seg(readStart=51, refIndex=1, refStart=1000, length=30, positive=true)
    val bp       = Breakpoint(segment1, segment2)
    (bp.leftPos, bp.leftPositive)   shouldBe (100, false)
    (bp.rightPos, bp.rightPositive) shouldBe (1000, true)

    // Same breakpoint observed with a "reverse" read where position in read is anti-correlated to position in genome
    val revSegment1 = seg(readStart=1, refIndex=1, refStart=1000, length=30, positive=false)
    val revSegment2 = seg(readStart=31,  refIndex=1, refStart=100,  length=50, positive=true)
    val rbp         = Breakpoint(revSegment1, revSegment2)
    (rbp.leftPos, rbp.leftPositive)   shouldBe (100, false)
    (rbp.rightPos, rbp.rightPositive) shouldBe (1000, true)
  }

  "Breakpoint.isCanonical" should "treat all strand combinations at the same position as canonical except (F, F)" in {
    val base = Breakpoint(leftRefIndex=0, leftPos=100, leftPositive=true, rightRefIndex=0, rightPos=100, rightPositive=true)

    // (T, T) at same position: canonical
    base.copy(leftPositive=true,  rightPositive=true).isCanonical  shouldBe true
    // (T, F) at same position: canonical (fixed-point of reversed)
    base.copy(leftPositive=true,  rightPositive=false).isCanonical shouldBe true
    // (F, T) at same position: canonical (fixed-point of reversed)
    base.copy(leftPositive=false, rightPositive=true).isCanonical  shouldBe true
    // (F, F) at same position: NOT canonical, reversed gives (T, T)
    base.copy(leftPositive=false, rightPositive=false).isCanonical shouldBe false
    base.copy(leftPositive=false, rightPositive=false).reversed.isCanonical shouldBe true
  }

  "Breakpoint.apply(segment, segment)" should "represent a tandem duplication sanely" in {
    // We want the same representation regardless of which piece of the duplication was sequenced more so
    // we test where the second segment is shorter, the same length as, and longer than the first
    Seq(-1, 0, 1).foreach { addend =>
      val segment1Length = 50
      val segment2Length = segment1Length + addend
      val segment1 = seg(readStart=1,  refIndex=1, refStart=100, length=segment1Length, positive=true)
      val segment2 = seg(readStart=51, refIndex=1, refStart=100, length=segment2Length, positive=true)
      val bp       = Breakpoint(segment1, segment2)
      (bp.leftPos,  bp.leftPositive)  shouldBe (100, false)
      (bp.rightPos, bp.rightPositive) shouldBe (149, false)
    }
  }
}