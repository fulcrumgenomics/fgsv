package com.fulcrumgenomics.sv

import com.fulcrumgenomics.util.Metric

/**
 * Represents a pileup of evidence (reads, read-pairs) for a breakpoint.  If `split_reads` is greater than
 * zero then the breakpoint location is supported at the sequence level and should be considered precise
 * though in the presence of repetitive sequence it may not be 100% accurate.  If `split_reads` is zero then
 * the only information comes from read-pairs and the breakpoint information should be considered imprecise.
 *
 * @param id            an ID assigned to the breakpoint that can be used to lookup supporting reads in the BAM.
 * @param left_contig   the contig of chromosome on which the left hand side of the breakpoint exists.
 * @param left_pos      the position (possibly imprecise) of the left-hand breakend.
 * @param left_strand   the strand of the left-hand breakend; sequence reads would traverse this strand
 *                      in order to arrive at the breakend and transit into the right-hand side of the breakpoint.
 * @param right_contig  the contig of chromosome on which the left hand side of the breakpoint exists.
 * @param right_pos     the position (possibly imprecise) of the right-hand breakend.
 * @param right_strand  the strand of the right-hand breakend;. sequence reads would continue reading onto
 *                      this strand after transiting the breakpoint from the left breakend
 * @param split_reads   the number of templates/inserts with split-read alignments that identified this breakpoint.
 * @param read_pairs    the number of templates/inserts with read-pair alignments (and without split-read alignments)
 *                      that identified this breakpoint.
 * @param total         the total number of templates/inserts that identified this breakpoint
 */
case class BreakpointPileup(id: String,
                            left_contig: String,
                            left_pos: Int,
                            left_strand: Char,
                            right_contig: String,
                            right_pos: Int,
                            right_strand: Char,
                            split_reads: Int,
                            read_pairs: Int,
                            total: Int
                           ) extends Metric

