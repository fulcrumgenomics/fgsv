package com.fulcrumgenomics.fgsv.tools

import com.fulcrumgenomics.fgsv.cmdline.{ClpGroups, SVTool}
import com.fulcrumgenomics.sopt.clp

@clp(group=ClpGroups.All, description=
  """
    |Trivial example tool to make sure build system is working.
  """)
class ExampleTool() extends SVTool {

  override def execute(): Unit = {
    System.err.println("Hello World!")
  }
}