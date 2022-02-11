package com.fulcrumgenomics.sv.util

import com.fulcrumgenomics.sv.UnitSpec
import com.fulcrumgenomics.sv.util.EvidenceType.SplitReadOppositeStrand

class EvidenceTypeTest extends UnitSpec {
  "EvidenceType.snakeName" should "convert to the name to snake-case" in {
    SplitReadOppositeStrand.snakeName shouldBe "split_read_opposite_strand"
  }
}
