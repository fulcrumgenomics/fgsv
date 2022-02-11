package com.fulcrumgenomics.sv.util

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.sv.UnitSpec
import com.fulcrumgenomics.sv.util.SegmentOrigin._
import com.fulcrumgenomics.testing.SamBuilder
import com.fulcrumgenomics.testing.SamBuilder.Strand

class AlignedSegmentTest extends UnitSpec {

  private def readEnd(start: Int, cigar: String, strand: Strand = SamBuilder.Plus, firstOfPair: Boolean = true)(implicit builder: SamBuilder): SamRecord = {
    val read = builder.addFrag(start=start, cigar=cigar, strand=strand).value
    read.paired = true
    read.firstOfPair = firstOfPair
    read
  }

  "AlignedSegment" should "be built from a SamRecord" in {
    implicit val builder: SamBuilder = new SamBuilder()

    // All bases mapped
    AlignedSegment(readEnd(start=1, cigar="100M")) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=0, start=1, end=100)
    )
    AlignedSegment(readEnd(start=1, cigar="100M", firstOfPair=false)) shouldBe AlignedSegment(
      origin=ReadTwo,
      readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=0, start=1, end=100),
    )
    AlignedSegment(readEnd(start=1, cigar="100M", strand=SamBuilder.Minus)) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=1, readEnd=100, positiveStrand=false, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=0, start=1, end=100)
    )

    // Leading soft-clipping
    AlignedSegment(readEnd(start=11, cigar="10S90M")) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=11, readEnd=100, positiveStrand=true, cigar=Cigar("10S90M"),
      range=GenomicRange(refIndex=0, start=11, end=100)
    )
    AlignedSegment(readEnd(start=11, cigar="10S90M", strand=SamBuilder.Minus)) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=1, readEnd=90, positiveStrand=false, cigar=Cigar("10S90M"),
      range=GenomicRange(refIndex=0, start=11, end=100)
    )

    // trailing soft-clipping
    AlignedSegment(readEnd(start=1, cigar="90M10S")) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=1, readEnd=90, positiveStrand=true, cigar=Cigar("90M10S"),
      range=GenomicRange(refIndex=0, start=1, end=90)
    )
    AlignedSegment(readEnd(start=1, cigar="90M10S", strand=SamBuilder.Minus)) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=11, readEnd=100, positiveStrand=false, cigar=Cigar("90M10S"),
      range=GenomicRange(refIndex=0, start=1, end=90)
    )

    // lots of clipping
    AlignedSegment(readEnd(start=10, cigar="1H4S90M6S2H")) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=6, readEnd=95, positiveStrand=true, cigar=Cigar("1H4S90M6S2H"),
      range=GenomicRange(refIndex=0, start=10, end=99)
    )
    AlignedSegment(readEnd(start=10, cigar="1H4S90M6S2H", strand=SamBuilder.Minus)) shouldBe AlignedSegment(
      origin=ReadOne,
      readStart=9, readEnd=98, positiveStrand=false, cigar=Cigar("1H4S90M6S2H"),
      range=GenomicRange(refIndex=0, start=10, end=99)
    )
  }

  "AlignmentBloc.segmentsFrom" should "create segments from a primary and one supplementals" in {
    val dummyRange = GenomicRange(refIndex=0, start=1, end=100)

    val primary = AlignedSegment(
      origin=ReadOne, readStart=1, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=dummyRange
    )
    val supplemental = AlignedSegment(
      origin=ReadOne, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange
    )

    // Keep both
    AlignedSegment.segmentsFrom(
      primary=primary, supplementals=Iterator(supplemental), readLength=100, minUniqueBasesToAdd=50
    ) should contain theSameElementsInOrderAs IndexedSeq(primary, supplemental)

    // Keep only the primary, as the supplemental doesn't add enough new bases
    AlignedSegment.segmentsFrom(
      primary=primary, supplementals=Iterator(supplemental), readLength=100, minUniqueBasesToAdd=51
    ) should contain theSameElementsInOrderAs IndexedSeq(primary)

    // Keep one supplemental
    AlignedSegment.segmentsFrom(
      primary=primary, supplementals=Iterator(supplemental, supplemental), readLength=100, minUniqueBasesToAdd=50
    ) should contain theSameElementsInOrderAs IndexedSeq(primary, supplemental)

    // A more complicated test case
    // Will be added, as we have 10 new bases
    val s1 = AlignedSegment(origin=ReadOne, readStart=40, readEnd=60, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Wont not be added, as we have only 9 new bases
    val s2 = AlignedSegment(origin=ReadOne, readStart=50, readEnd=69, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Will be added, as we have 10 new bases
    val s3 = AlignedSegment(origin=ReadOne, readStart=50, readEnd=70, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Wont not be added, as we have only 9 new bases
    val s4 = AlignedSegment(origin=ReadOne, readStart=92, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Will be added, as we have 10 new bases
    val s5 = AlignedSegment(origin=ReadOne, readStart=91, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Will not be added, as we have zero new bases
    val s6 = AlignedSegment(origin=ReadOne, readStart=91, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    AlignedSegment.segmentsFrom(
      primary=primary, supplementals=Iterator(s1, s2, s3, s4, s5, s6), readLength=100, minUniqueBasesToAdd=10
    ) should contain theSameElementsInOrderAs IndexedSeq(primary, s1, s3, s5)
  }

  "AlignedSegment.mergeReadSegments" should " merge read segments" in {
    val r1 = GenomicRange(refIndex=0, start=1, end=50) // overlaps nothing
    val r2 = GenomicRange(refIndex=0, start=100, end=150) // overlaps r1
    val r3 = GenomicRange(refIndex=0, start=125, end=175) // overlaps r2
    val r4 = GenomicRange(refIndex=1, start=1, end=1000) // overlaps nothing

    // No overlaps with one segment on R1 and two on R2
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignedSegment(origin=ReadTwo, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b3 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=r4)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1), r2Segments=IndexedSeq(b2, b3), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3)
    }

    // No overlaps with one segment on R1 and two on R2 due to strand
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignedSegment(origin=ReadTwo, readStart=51, readEnd=100, positiveStrand=false, cigar=Cigar.empty, range=r1)
      val b3 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=50, positiveStrand=false, cigar=Cigar.empty, range=r1)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1), r2Segments=IndexedSeq(b2, b3), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3)
    }

    // No overlaps with two segments on R1 and one on R2
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b2 = AlignedSegment(origin=ReadOne, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b3 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1), r2Segments=IndexedSeq(b2, b3), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3)
    }

    // No overlaps with two segments on R1 and two segments on R2
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=51, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignedSegment(origin=ReadOne, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b3 = AlignedSegment(origin=ReadTwo, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b4 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=r4)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1, b2), r2Segments=IndexedSeq(b3, b4), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3, b4)
    }

    // Single overlap on the last segment of R1 and the first segment of R2
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=51, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignedSegment(origin=ReadOne, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r2)
      val b3 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b23 = b2.merge(b3)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1, b2), r2Segments=IndexedSeq(b3), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b23)
    }

    // Single overlap on the last segment of R1 and the first segment of R2
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=51, positiveStrand=true, cigar=Cigar.empty, range=r2)
      val b2 = AlignedSegment(origin=ReadTwo, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b3 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b12 = b1.merge(b2)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1), r2Segments=IndexedSeq(b2, b3), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b12, b3)
    }

    // Length-two segment overlap
    {
      val b1 = AlignedSegment(origin=ReadOne, readStart=1, readEnd=30, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignedSegment(origin=ReadOne, readStart=31, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b3 = AlignedSegment(origin=ReadOne, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r2)
      val b4 = AlignedSegment(origin=ReadTwo, readStart=51, readEnd=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b5 = AlignedSegment(origin=ReadTwo, readStart=31, readEnd=50, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b6 = AlignedSegment(origin=ReadTwo, readStart=1, readEnd=30, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b24 = b2.merge(b4)
      val b35 = b3.merge(b5)
      AlignedSegment.mergeAlignedSegments(
        r1Segments=IndexedSeq(b1, b2, b3), r2Segments=IndexedSeq(b4, b5, b6), numOverlappingSegments=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b24, b35, b6)
    }
  }

  "SegmentOrigin.isPairedWith" should "return true if two segment origins can be considered a read pair" in {
    // Simple R1/R2 and R2/R1
    ReadOne.isPairedWith(ReadTwo) shouldBe true
    ReadTwo.isPairedWith(ReadOne) shouldBe true

    // Same read for the template, so not paired
    ReadOne.isPairedWith(ReadOne) shouldBe false
    ReadTwo.isPairedWith(ReadTwo) shouldBe false

    // One is R1/R2, the other is Both, so call it paired
    ReadOne.isPairedWith(Both) shouldBe true
    ReadTwo.isPairedWith(Both) shouldBe true
    Both.isPairedWith(ReadOne) shouldBe true
    Both.isPairedWith(ReadTwo) shouldBe true

    // Both are "Both", so it is compatible with R1/R2 (or R2/R1)
    Both.isPairedWith(Both) shouldBe true
  }
}
