package com.fulcrumgenomics.sv

import com.fulcrumgenomics.FgBioDef.PathToBam
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.{Matchers, OptionValues}

import java.nio.file.{Files, Path}

/** Base class for unit tests. */
trait UnitSpec extends AnyFlatSpec with Matchers with OptionValues {
  // Turn down HTSJDK logging
  htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)

  /** Creates a new temp file for use in testing that will be deleted when the VM exits. */
  protected def makeTempFile(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    path.toFile.deleteOnExit()
    path
  }

  /** Reads all the records from a SAM or BAM file into an indexed seq. */
  protected def readBamRecs(bam: PathToBam): IndexedSeq[SamRecord] = SamSource(bam).toIndexedSeq
}