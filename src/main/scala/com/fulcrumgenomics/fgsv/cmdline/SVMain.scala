package com.fulcrumgenomics.fgsv.cmdline

import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.FgBioMain
import com.fulcrumgenomics.vcf.api.VcfWriter

/**
  * Main program that loads everything up and runs the appropriate sub-command
  */
object SVMain {
  /** The main method */
  def main(args: Array[String]): Unit = new SVMain().makeItSoAndExit(args)
}

class SVMain extends FgBioMain {

  // HACK
  SamSource.DefaultUseAsyncIo = true
  SamWriter.DefaultUseAsyncIo = true
  VcfWriter.DefaultUseAsyncIo = true

   /** The name of the toolkit, used in printing usage and status lines. */
  override def name: String = "fgsv"

  /** The packages we wish to include in our command line **/
  override protected def packageList: List[String] =
    List[String]("com.fulcrumgenomics.fgsv")
}
