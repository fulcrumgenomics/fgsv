package com.fulcrumgenomics.sv

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.commons.util.StringUtil
import enumeratum.EnumEntry

import scala.collection.immutable

/** An enumeration over types of breakpoint evidence. */
object EvidenceType extends FgBioEnum[EvidenceType] {
  /** Evidence of the breakpoint was observed within aligned segments of a single read with split alignments. */
  case object SplitRead extends EvidenceType

  /** Evidence of the breakpoint was only observed _between_ reads in a read pair. */
  case object ReadPair extends EvidenceType

  override val values: immutable.IndexedSeq[EvidenceType] = findValues
}

sealed trait EvidenceType extends EnumEntry with Ordered[EvidenceType] {
  /** Provides ordering of evidence based on the order defined. */
  override def compare(that: EvidenceType): Int = EvidenceType.indexOf(this) - EvidenceType.indexOf(that)
  /** Returns teh snake-case name of the evidence type. */
  def snakeName: String = StringUtil.camelToGnu(this.entryName).replace('-', '_')
}