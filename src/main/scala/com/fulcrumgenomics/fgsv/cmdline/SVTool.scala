package com.fulcrumgenomics.fgsv.cmdline

import com.fulcrumgenomics.cmdline.FgBioTool
import com.fulcrumgenomics.commons.util.LazyLogging

/** All tools should extend this. */
trait SVTool extends FgBioTool with LazyLogging
