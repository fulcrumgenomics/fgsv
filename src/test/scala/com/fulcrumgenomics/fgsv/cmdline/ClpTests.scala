package com.fulcrumgenomics.fgsv.cmdline

import com.fulcrumgenomics.commons.util.CaptureSystemStreams
import com.fulcrumgenomics.fgsv.UnitSpec

/** Some basic test for the CLP classes. */
class ClpTests extends UnitSpec with CaptureSystemStreams {

  "SVTool" should "should print hello world" in {
    val (output, _, _) = captureItAll { () =>
      new SVMain().makeItSo("ExampleTool".split(' ')) shouldBe 0
    }
    output shouldBe "Hello World!\n"
  }
}
