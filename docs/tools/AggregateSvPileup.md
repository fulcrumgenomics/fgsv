---
title: AggregateSvPileup
---

# AggregateSvPileup

## Overview
**Group:** All tools

Merges nearby pileups of reads supporting putative breakpoints.

Takes as input the file of pileups produced by `SvPileup`. That file contains a list of breakpoints, each
consisting of a chromosome, position and strand for each side of the breakpoint, as well as quantified read support
for the breakpoint.

This tool merges sets of breakpoints that have their left sides on the same chromosome, their right sides on the
same chromosome, the same left and right strands, and their left and right positions both within a length
threshold. The merging behavior is transitive. For example, two breakpoints that are farther apart than the length
threshold can be merged if there is a third breakpoint that is close to both.

`SvPileup` distinguishes between two types of evidence for breakpoint events: split-read evidence, where a
breakpoint occurs within a mapped read, and read-pair evidence, where a breakpoint occurs in the unsequenced
insert between two paired reads. Currently this tool treats both forms of evidence equally, despite the
inaccuracy of positions reported by `SvPileup` for read-pair evidence.

If a bam file is provided, each aggregated pileup is annotated with the allele frequency at its left and right
breakends. Allele frequency is defined as the total number of templates supporting constituent breakpoints divided
by the count of the union of templates that cross any constituent breakends. In particular, paired templates that
straddle a breakend are considered to cross the breakend.

If a bed file of target regions is provided, each aggregated pileup is annotated with whether its left and
right sides overlap a target region.

The output file is a tab-delimited table with one record per aggregated cluster of pileups. Aggregated
pileups are reported with the minimum and maximum (inclusive) coordinates of all pileups in the cluster, a
possible putative structural variant event type supported by the pileups, and the sum of read support from all
pileups in the cluster.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|FilePath|Input text file of pileups generated by SvPileup|Required|1||
|bam|b|PathToBam|Bam file for allele frequency analysis. Must be the same file that was used as input to SvPileup.|Optional|1||
|flank|f|Int|If bam file is provided: distance upstream and downstream of aggregated breakpoint regions to search for mapped templates that overlap breakends. These are the templates that will be partitioned into those supporting the breakpoint vs. reading through it for the allele frequency calculation. Recommended to use at least the max read pair inner distance used by SvPileup.|Optional|1|1000|
|min-breakpoint-support||Int|If bam file is provided: minimum total number of templates supporting an aggregated breakpoint to report allele frequency. Supports speed improvement by avoiding querying and iterating over huge read pileups that contain insufficient support for a breakpoint to be considered interesting.|Optional|1|10|
|min-frequency||Double|If bam file is provided: minimum allele frequency to report. Supports speed improvement by avoiding iterating over huge read pileups that contain insufficient support for a breakpoint to be considered interesting.|Optional|1|0.001|
|targets-bed|t|FilePath|Optional bed file of target regions|Optional|1||
|output|o|FilePath|Output file|Required|1||
|max-dist|d|Int|Distance threshold below which to cluster breakpoints|Optional|1|10|

