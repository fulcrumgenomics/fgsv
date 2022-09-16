# Introduction

The following sections will help you to get started.

* [Installation](02_Installation.md)
* [Contributing](03_Contributing.md)
* [Metric Descriptions](04_Metrics.md)
* [Tools Descriptions](05_Tools.md)

## Overview

`fgsv` contains tools for gathering evidence for structural variants
from aligned reads. The `SvPileup` tool searches for split read mappings
and read pairs that map across breakpoints, emitting verbose information
similar to other "piluep" tools for small variant detection, but in this 
case for structural variation detection.  The `AggregateSvPileup` attempts
to aggregate information across "nearby" pileups, which is useful as often
the genomic start and end of a breakpoint is not always precise.  The tools
aim to be as sensitive as possible to find these evidence, but do neither
perform structural variation calling nor genotyping.
