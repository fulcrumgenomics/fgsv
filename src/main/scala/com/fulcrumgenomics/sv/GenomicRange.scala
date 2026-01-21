package com.fulcrumgenomics.sv

import htsjdk.samtools.util.CoordMath

/** A genomic interval and useful methods.
 *
 * @param refIndex the 0-based reference sequence index (as defined in the BAM/SAM header)
 * @param start    the 1-based start position (inclusive)
 * @param end      the 1-based end position (inclusive)
 */
case class GenomicRange(refIndex: Int, start: Int, end: Int) extends Ordered[GenomicRange] {
  require(start <= end, s"Start ($start) must be <= end ($end)")

  /** Returns true if the other interval overlaps this one (same reference, overlapping coordinates). */
  def overlaps(other: GenomicRange): Boolean = {
    this.refIndex == other.refIndex && CoordMath.overlaps(this.start, this.end, other.start, other.end)
  }

  /** Returns the union of the two genomic intervals.  The intervals must overlap. */
  def union(other: GenomicRange): GenomicRange = {
    require(this.overlaps(other), s"Can't union non-overlapping ranges: ${this} and ${other}")
    this.copy(start = Math.min(this.start, other.start), end = Math.max(this.end, other.end))
  }

  /** Compares to the other interval, returning less than 0 if this interval is before the other, 0 if the same, and > 0
   * if this interval is after the other. Orders by reference index, then start position, then end position. */
  def compare(that: GenomicRange): Int = {
    if (this.refIndex != that.refIndex) this.refIndex.compare(that.refIndex)
    else if (this.start != that.start) this.start.compare(that.start)
    else this.end.compare(that.end)
  }

  /** The length of the genomic range in bp. */
  def length: Int = this.end - this.start + 1

  override def toString: String =  s"${refIndex}:${start}-${end}"
}
