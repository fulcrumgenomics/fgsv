package com.fulcrumgenomics.fgsv.util

/** Stores the genomic break ends of a putative breakpoint
 *
 * @param id the identifier for this breakpoint
 * @param leftRefIndex the reference index (0-based) of the left side of the breakpoint
 * @param leftPos the genomic position (1-based) of the left side of the breakpoint
 * @param rightRefIndex  the reference index (0-based) of the right side of the breakpoint
 * @param rightPos the genomic position (1-based) of the right side of the breakpoint
 * @param sameStrand true if the sides are on the same strand, false otherwise
 * @param evidence the type of evidence for this breakpoint
 */
case class PutativeBreakpoint
(id: Long = 0,
 leftRefIndex: Int,
 leftPos: Int,
 rightRefIndex: Int,
 rightPos: Int,
 sameStrand: Boolean,
 evidence: EvidenceType
) extends Ordered[PutativeBreakpoint] {
  /** Returns true the other breakpoint has the same break ends.  This compares all fields except the evidence type. */
  def sameBreakEnds(that: PutativeBreakpoint): Boolean = {
    this.leftRefIndex == that.leftRefIndex &&
      this.leftPos == that.leftPos &&
      this.rightRefIndex == that.rightRefIndex &&
      this.rightPos == that.rightPos &&
      this.sameStrand == that.sameStrand
  }

  /** Defines an ordering over breakpoints, first by genomic coordinates of the left then right end, then evidence type. */
  def compare(that: PutativeBreakpoint): Int = {
    var result: Int = 0
    if (result == 0) result = this.leftRefIndex - that.leftRefIndex
    if (result == 0) result = this.leftPos - that.leftPos
    if (result == 0) result = this.rightRefIndex - that.rightRefIndex
    if (result == 0) result = this.rightPos - that.rightPos
    if (result == 0) result = EvidenceType.indexOf(this.evidence) - EvidenceType.indexOf(that.evidence)
    result
  }
}

object PutativeBreakpoint {
  /** Builds a breakpoint from two alignment blocks, with the left-side having the lower genomic coordinate. */
  def apply(b1: AlignmentBlock, b2: AlignmentBlock, evidence: EvidenceType): PutativeBreakpoint = {
    val (first, second) = if (b1.range <= b2.range) (b1, b2) else (b2, b1)
    PutativeBreakpoint(
      leftRefIndex  = first.range.refIndex,
      leftPos       = first.range.end,
      rightRefIndex = second.range.refIndex,
      rightPos      = second.range.start,
      sameStrand    = first.positiveStrand == second.positiveStrand,
      evidence      = evidence
    )
  }
}