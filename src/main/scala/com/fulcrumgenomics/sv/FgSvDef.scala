package com.fulcrumgenomics.sv

import com.fulcrumgenomics.FgBioDef.FilePath
import htsjdk.samtools.util.OverlapDetector
import htsjdk.tribble.AbstractFeatureReader
import htsjdk.tribble.bed.{BEDCodec, BEDFeature}

object FgSvDef {
  /** Builds an [[OverlapDetector]] frpom a BED file */
  def overlapDetectorFrom(bed: FilePath): OverlapDetector[BEDFeature] = {
    val bedReader = AbstractFeatureReader.getFeatureReader (bed.toAbsolutePath.toString, new BEDCodec (), false)
    val bedFeatures: java.util.List[BEDFeature] = bedReader.iterator ().toList
    bedReader.close()
    OverlapDetector.create (bedFeatures)
  }
}
