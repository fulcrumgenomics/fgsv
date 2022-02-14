package com.fulcrumgenomics.sv

import com.fulcrumgenomics.sv.EvidenceType.SplitReadOppositeStrand

class EvidenceTypeTest extends UnitSpec {
  "EvidenceType.snakeName" should "convert to the name to snake-case" in {
    SplitReadOppositeStrand.snakeName shouldBe "split_read_opposite_strand"
  }
}
