package com.fulcrumgenomics.sv.cmdline

import com.fulcrumgenomics.commons.util.CaptureSystemStreams
import com.fulcrumgenomics.sv.UnitSpec

/** Some basic test for the CLP classes. */
class ClpTests extends UnitSpec with CaptureSystemStreams {

  "SvTool" should "should print hello world" in {
    val (output, _, _) = captureItAll { () =>
      new SvMain().makeItSo("ExampleTool".split(' ')) shouldBe 0
    }
    output.contains("Hello World!") shouldBe true
  }
}
