package com.fulcrumgenomics.sv.cmdline

import com.fulcrumgenomics.cmdline.FgBioTool
import com.fulcrumgenomics.commons.util.LazyLogging

/** All tools should extend this. */
trait SvTool extends FgBioTool with LazyLogging
