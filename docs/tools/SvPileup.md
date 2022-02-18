---
title: SvPileup
---

# SvPileup

## Overview
**Group:** All tools

Collates a pileup of putative structural variant supporting reads.

Two output files will be created:

1. `<output-prefix>.txt`: a tab-delimited file describing SV pileups, one line per breakpiont event.
2. `<output-prefix>.bam`: a SAM/BAM file containing reads that contain SV breakpoint evidence annotated with SAM
  tags.  The `ev` SAM tag lists the type of evidence found, while the `be` list the unique breakpoint identifier
  output in (1) above.

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

