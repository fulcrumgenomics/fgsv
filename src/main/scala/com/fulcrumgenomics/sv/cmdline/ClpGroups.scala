package com.fulcrumgenomics.sv.cmdline

import com.fulcrumgenomics.sopt.cmdline.ClpGroup

/** Groups for organizing command line programs for display. */
object ClpGroups {

  class _All extends ClpGroup {
    override val name: String = "All tools"
    override val description: String = "All tools."
  }

  final val All = classOf[_All]
}
