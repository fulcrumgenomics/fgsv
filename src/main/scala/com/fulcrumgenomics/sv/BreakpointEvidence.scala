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
    new BreakpointEvidence(
      breakpoint = Breakpoint(from, into),
      evidence   = evidence,
      from       = (if (from.positiveStrand) from.right else from.left).toSet,
      into       = (if (into.positiveStrand) into.left else into.right).toSet
    )
  }
}


/** Stores the genomic break ends of a putative breakpoint
 *
 * @param breakpoint the breakpoint for which we have evidence
 * @param evidence the type of evidence for this breakpoint
 * @param from the [[SamRecord]](s) that provided evidence of this break point going "from" the breakpoint
 * @param into the [[SamRecord]](s) that provided evidence of this break point going "into" the breakpoint
 */
case class BreakpointEvidence(breakpoint: Breakpoint,
                              evidence: EvidenceType,
                              from: Set[SamRecord] = Set.empty,
                              into: Set[SamRecord] = Set.empty)
