package com.fulcrumgenomics.sv

import com.fulcrumgenomics.bam.api.SamRecord

object BreakpointEvidence {
  /**
   * Builds a breakpoint evidence from two aligned segments and an evidence type.
   *
   * The two segments must be given in such a way that a read traversing the breakpoint would
   * read through `from` in the genomic orientation given by `from` and then into `into` in the
   * orientation given by `into`.  E.g. if the two segments originate from a single (split-read)
   * alignment then `from` should come from earlier in the read in sequencing order.  Since the
   * segments could be from different reads in a mate pair, this cannot be validated/required()
   * in this method, but violating this assumption will lead to invalid breakpoints.
   *
   * The returned breakpoint will be canonicalized such that the `left` side of the breakpoint will have the
   * lower (or equal to) position on the genome vs. the `right` side.
   */
  def apply(from: AlignedSegment,
            into: AlignedSegment,
            evidence: EvidenceType): BreakpointEvidence = {
    new BreakpointEvidence(breakpoint=Breakpoint(from, into), evidence=evidence, recs=(from.recs++into.recs).toSet)
  }
}


/** Stores the genomic break ends of a putative breakpoint
 *
 * @param breakpoint the breakpoint for which we have evidence
 * @param evidence the type of evidence for this breakpoint
 * @param recs the [[SamRecord]](s) that provided evidence of this break point
 */
case class BreakpointEvidence(breakpoint: Breakpoint,
                              evidence: EvidenceType,
                              recs: Set[SamRecord] = Set.empty) extends Ordered[BreakpointEvidence] {

  /** Defines an ordering over breakpoints, first by genomic coordinates of the left then right end, then evidence type.
   *  Ignores the origin of the breakpoint and evidence */
  def compare(that: BreakpointEvidence): Int = {
    var result              = this.breakpoint.compare(that.breakpoint)
    if (result == 0) result = this.evidence.compare(that.evidence)
    result
  }
}
