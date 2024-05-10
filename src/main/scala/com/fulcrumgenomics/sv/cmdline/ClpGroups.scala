package com.fulcrumgenomics.sv.cmdline

import com.fulcrumgenomics.sopt.cmdline.ClpGroup

/** Groups for organizing command line programs for display. */
object ClpGroups {

  class _BreakpointAndSv extends ClpGroup {
    override val name: String = "Breakpoint and SV Tools"
    override val description: String = "Primary tools for calling and transforming breakpoints and SVs."
  }

  class _Utilities  extends ClpGroup {
    override val name: String = "Utility Tools"
    override val description: String = "Helper tools for working with breakpoint or SV data."
  }

  final val BreakpointAndSv = classOf[_BreakpointAndSv]
  final val Utilities = classOf[_Utilities]
}
