package com.fulcrumgenomics.sv

import org.scalatest.{FlatSpec, Matchers, OptionValues}

import java.nio.file.{Files, Path}

/** Base class for unit tests. */
class UnitSpec extends FlatSpec with Matchers with OptionValues {

  /** Creates a new temp file for use in testing that will be deleted when the VM exits. */
  protected def makeTempFile(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    path.toFile.deleteOnExit()
    path
  }

}
