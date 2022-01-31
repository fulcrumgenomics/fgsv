package com.fulcrumgenomics.fgsv.util

import com.fulcrumgenomics.fgsv.UnitSpec
import com.fulcrumgenomics.fgsv.util.EvidenceType.SplitReadOppositeStrand

class EvidenceTypeTest extends UnitSpec {
  "EvidenceType.snakeName" should "convert to the name to snake-case" in {
    SplitReadOppositeStrand.snakeName shouldBe "split_read_opposite_strand"
  }
}
