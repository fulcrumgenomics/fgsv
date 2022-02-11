package com.fulcrumgenomics.sv

import htsjdk.samtools.util.CoordMath

/** A genomic interval and useful methods. */
case class GenomicRange(refIndex: Int, start: Int, end: Int) extends Ordered[GenomicRange] {
  require(start <= end)

  /** Returns if the other interval overlaps this one. */
  def overlaps(other: GenomicRange): Boolean = {
    this.refIndex == other.refIndex && CoordMath.overlaps(this.start, this.end, other.start, other.end)
  }

  /** Returns the union of the two genomic intervals.  The intervals must overlap. */
  def union(other: GenomicRange): GenomicRange = {
    require(this.overlaps(other))
    this.copy(start = Math.min(this.start, other.start), end = Math.max(this.end, other.end))
  }

  /** Compares to the other interval, return less than 0 if this interval is before the other, 0 if the same, and > 0
   * if this interval is after the other. */
  def compare(that: GenomicRange): Int = {
    if (this.refIndex != that.refIndex) this.refIndex - that.refIndex
    else if (this.start != that.start) this.start - that.start
    else this.end - that.end
  }

  /** The length of the genomic range in bp. */
  def length: Int = this.end - this.start + 1
}
