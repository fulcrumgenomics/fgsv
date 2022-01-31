package com.fulcrumgenomics.sv.tools

import com.fulcrumgenomics.sv.cmdline.{ClpGroups, SvTool}
import com.fulcrumgenomics.sopt.clp

@clp(group=ClpGroups.All, description=
  """
    |Trivial example tool to make sure build system is working.
  """)
class ExampleTool() extends SvTool {

  override def execute(): Unit = {
    System.err.println("Hello World!")
  }
}