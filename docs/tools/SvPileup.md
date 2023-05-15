---
title: SvPileup
---

# SvPileup

## Overview
**Group:** All tools

Collates a pileup of putative structural variant supporting reads.

## Outputs

Two output files will be created:

1. `<output-prefix>.txt`: a tab-delimited file describing SV pileups, one line per breakpoint event.  The returned
   breakpoint will be canonicalized such that the "left" side of the breakpoint will have the lower (or equal to)
   position on the genome vs. the "right"s side.
2. `<output-prefix>.bam`: a SAM/BAM file containing reads that contain SV breakpoint evidence annotated with SAM
  tag.

The `be` SAM tag contains a comma-delimited list of breakpoints to which a given read belongs.  Each element is
a semi-colon delimited, with four fields:

1. The unique breakpoint identifier (same identifier found in the tab-delimited output).
2. Either "left" or "right, corresponding to if the read shows evidence of the genomic left or right side of the
   breakpoint as found in the breakpoint file (i.e. `left_pos` or `right_pos`).
3. Either "from" or "into", such that when traversing the breakpoint would read through "from" and then into
   "into" in the sequencing order of the read pair.  For a split-read alignment, the "from" contains the aligned
   portion of the read that comes from earlier in the read in sequencing order.  For an alignment of a read-pair
   spanning the breakpoint, then "from" should be read-one of the pair and "into" should be read-two of the pair.
4. The type of breakpoint evidence: either "split_read" for observations of an aligned segment of a single read
   with split alignments, or "read_pair" for observations _between_ reads in a read pair.

## Example output

The following shows two breakpoints:

```
id left_contig left_pos left_strand right_contig right_pos right_strand split_reads read_pairs total
 1        chr1      100           +         chr2       200            -           1          0     1
 2        chr2      150           -         chr3       500            +           1          0     1
```

Consider a single fragment read that maps across both the above two breakpoints, so has three split-read
alignments.  The first alignment maps on the left side of breakpoint #1, the second alignment maps to both the
right side of breakpoint #1 and the left-side of breakpoint #2, and the third alignment maps to the right side of
breakpoint #2. The SAM records would be as follows:

```
r1    0 chr1  50 60   50M100S ... be:Z:1;left;from;split_read
r1 2064 chr2 150 60 50S50M50S ... be:Z:1;right;into;split_read,2;left;from;split_read
r1 2048 chr3 500 60   100S50M ... be:Z:2;right;into;split_read
```

## Algorithm Overview

Putative breakpoints are identified by examining the alignments for each template. The alignments are transformed
into aligned segments in the order they were sequenced.  Each aligned segment represents the full genomic span of
the mapped bases.  This is performed first for the primary alignments.  Next, supplementary alignments are added
only if they map read bases that have not been previously covered by other alignments (see
`--min-unique-bases-to-add`).  This is iteratively performed until supplementary alignments have been exhausted.

Next, aligned segments that have overlapping genomic mapped bases are merged into a single aligned
segment.  In this case, the two or more read mappings merged are associated with either the left side or right
side of that aligned segment, controlled by examining how close to the end of the new aligned segment the given
read mapping occurs (see `--slop` option).  This used to identify which reads traverse "from" and "into" the
breakpoint as described above.

Finally, pairs of adjacent aligned segments are examined for evidence of a breakpoint, such that genomic distance
between them beyond either `--max-read-pair-inner-distance` for aligned segments from different read pairs, or
`--max-aligned-segment-inner-distance` for aligned segments from the same read in a pair (i.e. split-read mapping).
Split read evidence will be returned in favor of across-read-pair evidence when both are present.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|The input query sorted or grouped BAM|Required|1||
|output|o|PathPrefix|The output path prefix|Required|1||
|max-read-pair-inner-distance|d|Int|The maximum _inner_ distance for normal read pair|Optional|1|1000|
|max-aligned-segment-inner-distance|D|Int|The maximum _inner_ distance between two segments of a split read mapping|Optional|1|100|
|min-primary-mapping-quality|q|Int|The minimum mapping quality for primary alignments|Optional|1|30|
|min-supplementary-mapping-quality|Q|Int|The minimum mapping quality for supplementary alignments|Optional|1|18|
|min-unique-bases-to-add|b|Int|The minimum # of uncovered query bases needed to add a supplemental alignment|Optional|1|20|
|slop|s|Int|The number of bases of slop to allow when determining which records to track for the left or right side of an aligned segment when merging segments.|Optional|1|5|
|targets-bed|t|FilePath|Optional bed file of target regions|Optional|1||
|targets-bed-requirement|T|Requirement|Requirement on if each side of the breakpoint must overlap a target.  Will always annotate each side of the breakpoint.|Optional|1|AnnotateOnly|

