package com.fulcrumgenomics.fgsv.util

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.commons.util.StringUtil
import enumeratum.EnumEntry

import scala.collection.immutable

/** An enumeration over types of SV evidence. */
object EvidenceType extends FgBioEnum[EvidenceType] {

  /** Evidence from one end of a read mapping to a different contig than the other. */
  case object ReadPairInterContig extends EvidenceType
  /** Evidence from one end of a read mapping to a the same contig as the other, but too far apart. */
  case object ReadPairIntraContig extends EvidenceType
  /** Evidence from a reverse-forward read pair orientation. */
  case object ReadPairReverseForward extends EvidenceType
  /** Evidence from a tandem read pair orientation. */
  case object ReadPairTandem extends EvidenceType
  /** Evidence from a split read mapping where adjacent blocks map to different contigs. */
  case object SplitReadInterContig extends EvidenceType
  /** Evidence from a split read mapping where adjacent blocks map to the same contig, but too far apart. */
  case object SplitReadIntraContig extends EvidenceType
  /** Evidence from a split read mapping where adjacent blocks map to opposite genomic strands. */
  case object SplitReadOppositeStrand extends EvidenceType

  override val values: immutable.IndexedSeq[EvidenceType] = findValues
}

sealed trait EvidenceType extends EnumEntry with Ordered[EvidenceType] {
  /** Provides ordering of evidence based on the order defined. */
  override def compare(that: EvidenceType): Int = EvidenceType.indexOf(this) - EvidenceType.indexOf(that)
  /** Returns teh snake-case name of the evidence type. */
  def snakeName: String = StringUtil.camelToGnu(this.entryName).replace('-', '_')
}