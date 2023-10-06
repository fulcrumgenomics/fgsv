---
title: FilterAndMerge
---

# FilterAndMerge

## Overview
**Group:** All tools

Filters and merges SVPileup output.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|FilePath|The input pileup file from SvPileup|Required|1||
|output|o|FilePath|The output filtered and merged SvPileup file|Required|1||
|dict|d|PathToSequenceDictionary|The path to the reference sequence dictionary.|Required|1||
|min-pre|m|Int|The minimum # of observations to examine an input site|Optional|1|1|
|min-post|M|Int|The minimum # of observations to output a site|Optional|1|1|
|slop|s|Int|The maximum # bases between a breakend across adjacent sites|Optional|1|0|

