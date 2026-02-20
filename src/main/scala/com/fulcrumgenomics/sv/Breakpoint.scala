package com.fulcrumgenomics.sv

import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.samtools.util.Interval

object Breakpoint {
  /**
   * Builds a breakpoint from two aligned segments.
   *
   * The two segments must be given in such a way that a read traversing the breakpoint would
   * read through `from` in the genomic orientation given by `from` and then into `into` in the
   * orientation given by `into`.  E.g. if the two segments originate from a single (split-read)
   * alignment then `from` should come from earlier in the read in sequencing order.  Since the
   * segments could be from different reads in a mate pair, this cannot be validated/required()
   * in this method, but violating this assumption will lead to invalid breakpoints.
   *
   * @param from the segment earlier in sequencing order
   * @param into the segment later in sequencing order
   * @param canonicalize canonicalize the breakpoint such that the `left` side of the breakpoint will have the lower
   *                     (or equal to) position on the genome vs. the `right` side.
   */
  def apply(from: AlignedSegment, into: AlignedSegment, canonicalize: Boolean = true): Breakpoint = {
    val bp = Breakpoint(
      leftRefIndex  = from.range.refIndex,
      leftPos       = if (from.positiveStrand) from.range.end else from.range.start,
      leftPositive  = from.positiveStrand,
      rightRefIndex = into.range.refIndex,
      rightPos      = if (into.positiveStrand) into.range.start else into.range.end,
      rightPositive = into.positiveStrand
    )

    if (!canonicalize || bp.isCanonical) bp else bp.reversed
  }

  /**
   * An alternative ordering for breakpoints that orders by left/right chrom, then left/right pos etc.
   * so that breakpoints between approximately the same positions on the same chromosomes group together.
   */
  val PairedOrdering: Ordering[Breakpoint] = new Ordering[Breakpoint] {
    override def compare(x: Breakpoint, y: Breakpoint): Int = {
      var result              = x.leftRefIndex.compare(y.leftRefIndex)
      if (result == 0) result = x.rightRefIndex.compare(y.rightRefIndex)
      if (result == 0) result = x.leftPos.compare(y.leftPos)
      if (result == 0) result = x.rightPos.compare(y.rightPos)
      if (result == 0) result = x.leftPositive.compare(y.leftPositive)
      if (result == 0) result = x.rightPositive.compare(y.rightPositive)
      result
    }
  }
}


/** Stores the genomic break ends of a possible breakpoint
 *
 * @param leftRefIndex the reference index (0-based) of the left side of the breakpoint
 * @param leftPos the genomic position (1-based) of the left side of the breakpoint
 * @param leftPositive true if the sequence on the left hand side of the breakpoint is on the positive genomic strand
 * @param rightRefIndex  the reference index (0-based) of the right side of the breakpoint
 * @param rightPos the genomic position (1-based) of the right side of the breakpoint
 * @param rightPositive true if the sequence on the right hand side of the breakpoint is on the positive genomic strand
 */
case class Breakpoint(leftRefIndex: Int,
                      leftPos: Int,
                      leftPositive: Boolean,
                      rightRefIndex: Int,
                      rightPos: Int,
                      rightPositive: Boolean,
                     ) extends Ordered[Breakpoint] {

  @inline final def leftNegative: Boolean = !leftPositive
  @inline final def rightNegative: Boolean = !rightPositive

  /** Orders breakpoints by their left genomic position then their right genomic position. */
  override def compare(that: Breakpoint): Int = {
    var result              = this.leftRefIndex.compare(that.leftRefIndex)
    if (result == 0) result = this.leftPos.compare(that.leftPos)
    if (result == 0) result = this.leftPositive.compare(that.leftPositive)
    if (result == 0) result = this.rightRefIndex.compare(that.rightRefIndex)
    if (result == 0) result = this.rightPos.compare(that.rightPos)
    if (result == 0) result = this.rightPositive.compare(that.rightPositive)
    result
  }

  /**
   * Returns an inverted copy of the breakpoint where the left and right sides of the breakpoint have been
   * swapped and the strand information has been flipped.  E.g. if the breakpoint were to specify
   * `chr1:500:F > chr1:100:F`, inverting it would yield `chr1:100:R > chr1:500R` to maintain consistency.
   */
  def reversed: Breakpoint = copy(
    leftRefIndex  = rightRefIndex,
    leftPos       = rightPos,
    leftPositive  = !rightPositive,
    rightRefIndex = leftRefIndex,
    rightPos      = leftPos,
    rightPositive = !leftPositive
  )

  /** Returns true if the representation of the breakpoint is canonical, with the left hand side of the break
   * earlier on the genome than the right hand side. */
  def isCanonical: Boolean = (leftRefIndex < rightRefIndex) ||
    (leftRefIndex == rightRefIndex && leftPos < rightPos) ||
    (leftRefIndex == rightRefIndex && leftPos == rightPos && leftPositive)

  /** Returns an [[Interval]] for the left side of this breakpoint */
  def leftInterval(dict: SequenceDictionary): Interval = new Interval(
    dict(leftRefIndex).name,
    leftPos,
    leftPos
  )

  /** Returns an [[Interval]] for the right side of this breakpoint */
  def rightInterval(dict: SequenceDictionary): Interval = new Interval(
    dict(rightRefIndex).name,
    rightPos,
    rightPos
  )

}
