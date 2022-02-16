package com.fulcrumgenomics.sv

import com.fulcrumgenomics.sv.EvidenceType.SplitRead

class EvidenceTypeTest extends UnitSpec {
  "EvidenceType.snakeName" should "convert to the name to snake-case" in {
    SplitRead.snakeName shouldBe "split_read"
  }
}
