package com.fulcrumgenomics.fgsv.util

import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.fgsv.UnitSpec
import com.fulcrumgenomics.fgsv.util.BlockOrigin._
import com.fulcrumgenomics.testing.SamBuilder
import com.fulcrumgenomics.testing.SamBuilder.Strand

class AlignmentBlockTest extends UnitSpec {

  private def readEnd(start: Int, cigar: String, strand: Strand = SamBuilder.Plus, firstOfPair: Boolean = true)(implicit builder: SamBuilder): SamRecord = {
    val read = builder.addFrag(start=start, cigar=cigar, strand=strand).value
    read.paired = true
    read.firstOfPair = firstOfPair
    read
  }

  "AlignmentBlock" should "be built from a SamRecord" in {
    implicit val builder: SamBuilder = new SamBuilder()

    // All bases mapped
    AlignmentBlock(readEnd(start=1, cigar="100M")) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=1, end=100, positiveStrand=true, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=0, start=1, end=100)
    )
    AlignmentBlock(readEnd(start=1, cigar="100M", firstOfPair=false)) shouldBe AlignmentBlock(
      origin=ReadTwo,
      start=1, end=100, positiveStrand=true, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=0, start=1, end=100),
    )
    AlignmentBlock(readEnd(start=1, cigar="100M", strand=SamBuilder.Minus)) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=1, end=100, positiveStrand=false, cigar=Cigar("100M"),
      range=GenomicRange(refIndex=0, start=1, end=100)
    )

    // Leading soft-clipping
    AlignmentBlock(readEnd(start=11, cigar="10S90M")) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=11, end=100, positiveStrand=true, cigar=Cigar("10S90M"),
      range=GenomicRange(refIndex=0, start=11, end=100)
    )
    AlignmentBlock(readEnd(start=11, cigar="10S90M", strand=SamBuilder.Minus)) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=1, end=90, positiveStrand=false, cigar=Cigar("10S90M"),
      range=GenomicRange(refIndex=0, start=11, end=100)
    )

    // trailing soft-clipping
    AlignmentBlock(readEnd(start=1, cigar="90M10S")) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=1, end=90, positiveStrand=true, cigar=Cigar("90M10S"),
      range=GenomicRange(refIndex=0, start=1, end=90)
    )
    AlignmentBlock(readEnd(start=1, cigar="90M10S", strand=SamBuilder.Minus)) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=11, end=100, positiveStrand=false, cigar=Cigar("90M10S"),
      range=GenomicRange(refIndex=0, start=1, end=90)
    )

    // lots of clipping
    AlignmentBlock(readEnd(start=10, cigar="1H4S90M6S2H")) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=6, end=95, positiveStrand=true, cigar=Cigar("1H4S90M6S2H"),
      range=GenomicRange(refIndex=0, start=10, end=99)
    )
    AlignmentBlock(readEnd(start=10, cigar="1H4S90M6S2H", strand=SamBuilder.Minus)) shouldBe AlignmentBlock(
      origin=ReadOne,
      start=9, end=98, positiveStrand=false, cigar=Cigar("1H4S90M6S2H"),
      range=GenomicRange(refIndex=0, start=10, end=99)
    )
  }

  "AlignmentBloc.blocksFrom" should "create blocks from a primary and one supplementals" in {
    val dummyRange = GenomicRange(refIndex=0, start=1, end=100)

    val primary = AlignmentBlock(
      origin=ReadOne, start=1, end=50, positiveStrand=true, cigar=Cigar.empty, range=dummyRange
    )
    val supplemental = AlignmentBlock(
      origin=ReadOne, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange
    )

    // Keep both
    AlignmentBlock.blocksFrom(
      primary=primary, supplementals=Iterator(supplemental), readLength=100, minUniqueBasesToAdd=50
    ) should contain theSameElementsInOrderAs IndexedSeq(primary, supplemental)

    // Keep only the primary, as the supplemental doesn't add enough new bases
    AlignmentBlock.blocksFrom(
      primary=primary, supplementals=Iterator(supplemental), readLength=100, minUniqueBasesToAdd=51
    ) should contain theSameElementsInOrderAs IndexedSeq(primary)

    // Keep one supplemental
    AlignmentBlock.blocksFrom(
      primary=primary, supplementals=Iterator(supplemental, supplemental), readLength=100, minUniqueBasesToAdd=50
    ) should contain theSameElementsInOrderAs IndexedSeq(primary, supplemental)

    // A more complicated test case
    // Will be added, as we have 10 new bases
    val s1 = AlignmentBlock(origin=ReadOne, start=40, end=60, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Wont not be added, as we have only 9 new bases
    val s2 = AlignmentBlock(origin=ReadOne, start=50, end=69, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Will be added, as we have 10 new bases
    val s3 = AlignmentBlock(origin=ReadOne, start=50, end=70, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Wont not be added, as we have only 9 new bases
    val s4 = AlignmentBlock(origin=ReadOne, start=92, end=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Will be added, as we have 10 new bases
    val s5 = AlignmentBlock(origin=ReadOne, start=91, end=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    // Will not be added, as we have zero new bases
    val s6 = AlignmentBlock(origin=ReadOne, start=91, end=100, positiveStrand=true, cigar=Cigar.empty, range=dummyRange)
    AlignmentBlock.blocksFrom(
      primary=primary, supplementals=Iterator(s1, s2, s3, s4, s5, s6), readLength=100, minUniqueBasesToAdd=10
    ) should contain theSameElementsInOrderAs IndexedSeq(primary, s1, s3, s5)
  }

  "AlignmentBlock.mergeReadBlocks" should " merge read blocks" in {
    val r1 = GenomicRange(refIndex=0, start=1, end=50) // overlaps nothing
    val r2 = GenomicRange(refIndex=0, start=100, end=150) // overlaps r1
    val r3 = GenomicRange(refIndex=0, start=125, end=175) // overlaps r2
    val r4 = GenomicRange(refIndex=1, start=1, end=1000) // overlaps nothing

    // No overlaps with one block on R1 and two on R2
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignmentBlock(origin=ReadTwo, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b3 = AlignmentBlock(origin=ReadTwo, start=1, end=50, positiveStrand=true, cigar=Cigar.empty, range=r4)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1), r2Blocks=IndexedSeq(b2, b3), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3)
    }

    // No overlaps with one block on R1 and two on R2 due to strand
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignmentBlock(origin=ReadTwo, start=51, end=100, positiveStrand=false, cigar=Cigar.empty, range=r1)
      val b3 = AlignmentBlock(origin=ReadTwo, start=1, end=50, positiveStrand=false, cigar=Cigar.empty, range=r1)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1), r2Blocks=IndexedSeq(b2, b3), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3)
    }

    // No overlaps with two blocks on R1 and one on R2
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=50, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b2 = AlignmentBlock(origin=ReadOne, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b3 = AlignmentBlock(origin=ReadTwo, start=1, end=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1), r2Blocks=IndexedSeq(b2, b3), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3)
    }

    // No overlaps with two blocks on R1 and two blocks on R2
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=51, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignmentBlock(origin=ReadOne, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b3 = AlignmentBlock(origin=ReadTwo, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b4 = AlignmentBlock(origin=ReadTwo, start=1, end=50, positiveStrand=true, cigar=Cigar.empty, range=r4)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1, b2), r2Blocks=IndexedSeq(b3, b4), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b2, b3, b4)
    }

    // Single overlap on the last block of R1 and the first block of R2
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=51, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignmentBlock(origin=ReadOne, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r2)
      val b3 = AlignmentBlock(origin=ReadTwo, start=1, end=100, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b23 = b2.merge(b3)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1, b2), r2Blocks=IndexedSeq(b3), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b23)
    }

    // Single overlap on the last block of R1 and the first block of R2
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=51, positiveStrand=true, cigar=Cigar.empty, range=r2)
      val b2 = AlignmentBlock(origin=ReadTwo, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b3 = AlignmentBlock(origin=ReadTwo, start=1, end=50, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b12 = b1.merge(b2)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1), r2Blocks=IndexedSeq(b2, b3), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b12, b3)
    }

    // Length-two block overlap
    {
      val b1 = AlignmentBlock(origin=ReadOne, start=1, end=30, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b2 = AlignmentBlock(origin=ReadOne, start=31, end=50, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b3 = AlignmentBlock(origin=ReadOne, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r2)
      val b4 = AlignmentBlock(origin=ReadTwo, start=51, end=100, positiveStrand=true, cigar=Cigar.empty, range=r1)
      val b5 = AlignmentBlock(origin=ReadTwo, start=31, end=50, positiveStrand=true, cigar=Cigar.empty, range=r3)
      val b6 = AlignmentBlock(origin=ReadTwo, start=1, end=30, positiveStrand=true, cigar=Cigar.empty, range=r4)
      val b24 = b2.merge(b4)
      val b35 = b3.merge(b5)
      AlignmentBlock.mergeReadBlocks(
        r1Blocks=IndexedSeq(b1, b2, b3), r2Blocks=IndexedSeq(b4, b5, b6), numOverlappingBlocks=1
      ) should contain theSameElementsInOrderAs IndexedSeq(b1, b24, b35, b6)
    }
  }

  "BlockOrigin.isPairedWith" should "return true if two block origins can be considered a read pair" in {
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
