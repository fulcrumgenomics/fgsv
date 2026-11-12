---
title: RescueSvJunction
---

# RescueSvJunction

## Overview
**Group:** Breakpoint and SV Tools

Experimental soft-clip rescue for one supplied, oriented junction. Compares a read against
the synthetic junction and both local reference alternatives, requiring aligned anchors on
both sides. Reports qualifying realigned support, including existing discovery supporters;
consumers must union template identities with original evidence, not add counts. Does not
modify the BAM, emit SA tags, discover junctions, or update pileup counts or frequencies.

Eligible reads must be mapped primary alignments, not secondary, supplementary, duplicate,
or vendor-QC-failed, and meet --min-mapq. A terminal soft clip must face the supplied endpoint
within --slop. Both the observed clip and the immediately adjacent contiguous aligned source
segment must meet --min-anchor; the final junction alignment must also meet --min-anchor on
each side. Hard clipping does not count as observed support. Equivalent adjacent CIGAR
operations are coalesced for eligibility, and padding does not interrupt a source anchor.
SAM SEQ '=' is expanded from the original aligned reference positions; '=' in soft clips
or insertions has no reference base and must be replaced with literal bases before rescue.

Output is a headered metrics table, including when no templates qualify. Each row identifies
the supplied endpoints and the actual synthetic-junction boundary after genomic-edge clamping.
Coordinates are 1-based inclusive. Query coordinates use the orientation in query_reversed;
target coordinates refer to the oriented synthetic junction. Endpoint strand conventions
match SvPileup's left/from and right/into endpoints. For a negative strand use --left-strand=-
or --right-strand=-. Deduplication applies only within this supplied junction.
For example, with flank 200, --left-position 600 --left-strand=+ --right-position 401
--right-strand=+ joins reference bases 401..600 to 401..600 at target boundary 200.
With --left-strand=- --right-strand=- at the same endpoints, it joins reverse-complemented
600..799 to reverse-complemented 202..401, also at boundary 200. Contigs are supplied separately.

The representative per template is selected by descending score advantage, junction score,
and minimum anchor, then ascending read number, genomic/alignment coordinates, CIGARs and
query orientation (false first), with remaining diagnostics breaking ties. This does not
override rejection of ambiguous optimal alignment placements.

Smith-Waterman scores are fixed: canonical match +2, mismatch/noncanonical base -3,
gap open -6 and gap extension -1. Anchors count canonical aligned columns, including
mismatches; gaps and ambiguous bases do not contribute. Scores compare only the supplied
local alternatives, not all genomic mappings, and are not quality-aware genotype likelihoods
or supplementary mapping qualities. The flank must cover each otherwise eligible read's
observed length plus slop; increase --flank for longer reads. Genomic edges are clamped,
not compensated by truncating reads.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Indexed coordinate-sorted original BAM|Required|1||
|reference|r|FilePath|Indexed reference FASTA used for the BAM|Required|1||
|output|o|FilePath|Output qualifying-template metrics; existing discovery supporters are included|Required|1||
|left-contig||String|Left/from reference contig|Required|1||
|left-position||Int|Left/from aligned endpoint, 1-based inclusive|Required|1||
|left-strand||String|Left/from strand, + or -|Required|1||
|right-contig||String|Right/into reference contig|Required|1||
|right-position||Int|Right/into aligned endpoint, 1-based inclusive|Required|1||
|right-strand||String|Right/into strand, + or -|Required|1||
|min-anchor||Int|Minimum aligned bases on each side of the junction|Optional|1|10|
|slop||Int|Maximum distance from a soft-clip boundary to its candidate endpoint|Optional|1|5|
|min-mapq||Int|Minimum primary alignment mapping quality|Optional|1|30|
|min-score-delta||Int|Minimum junction score advantage over both local reference alternatives|Optional|1|10|
|min-identity||Double|Minimum junction-alignment identity|Optional|1|0.9|
|min-query-coverage||Double|Minimum fraction of observed query bases included in the junction alignment|Optional|1|0.9|
|flank||Int|Reference flank and BAM query padding; must cover eligible read length plus slop|Optional|1|250|
