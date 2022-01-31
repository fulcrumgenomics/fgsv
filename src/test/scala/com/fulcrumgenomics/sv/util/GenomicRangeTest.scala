package com.fulcrumgenomics.sv.util

import com.fulcrumgenomics.FgBioDef.forloop
import com.fulcrumgenomics.sv.UnitSpec

class GenomicRangeTest extends UnitSpec {
  private val r1 = GenomicRange(refIndex=1, start=1, end=10)
  private val r2 = GenomicRange(refIndex=1, start=10, end=20)
  private val r3 = GenomicRange(refIndex=1, start=11, end=21)
  private val r4 = GenomicRange(refIndex=2, start=1, end=100)

  "GenomicRange.overlaps" should "return true if two ranges overlap" in {
    r1.overlaps(r2) shouldBe true
    r2.overlaps(r1) shouldBe true

    r1.overlaps(r3) shouldBe false
    r3.overlaps(r1) shouldBe false

    r1.overlaps(r4) shouldBe false
    r4.overlaps(r1) shouldBe false
  }

  "GenomicRange.union" should "union two overlapping ranges" in {
    r1.union(r2) shouldBe GenomicRange(refIndex=1, start=1, end=20)
    r2.union(r1) shouldBe GenomicRange(refIndex=1, start=1, end=20)

    r2.union(r3) shouldBe GenomicRange(refIndex=1, start=10, end=21)
    r3.union(r2) shouldBe GenomicRange(refIndex=1, start=10, end=21)
  }

  "Genomic.range" should "compare two ranges" in {
    val ranges = IndexedSeq(r1, r2, r3, r4)
    forloop(from=0, until=ranges.length) { leftIndex =>
      val left = ranges(leftIndex)
      left.compare(left) shouldBe 0
      forloop(from=leftIndex + 1, until=ranges.length) { rightIndex =>
        val right = ranges(rightIndex)
        left.compare(right) should be < 0
        right.compare(left) should be > 0
      }
    }
  }
}
