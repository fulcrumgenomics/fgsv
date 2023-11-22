# fgsv

[![Bioconda][bioconda-badge-link]][bioconda-link]
[![Build Status][github-badge]][github-link]
[![Language][scala-badge]][scala-link]
[![License][license-badge]][license-link]

[bioconda-badge-link]: https://img.shields.io/conda/dn/bioconda/fgsv.svg?label=Bioconda
[bioconda-link]:       http://bioconda.github.io/recipes/fgsv/README.html
[github-badge]:        https://github.com/fulcrumgenomics/fgsv/actions/workflows/unittests.yaml/badge.svg?branch=main
[github-link]:         https://github.com/fulcrumgenomics/fgsv/actions/workflows/unittests.yaml
[scala-badge]:         https://img.shields.io/badge/language-scala-c22d40.svg
[scala-link]:          https://www.scala-lang.org/
[license-badge]:       https://img.shields.io/badge/license-MIT-blue.svg
[license-link]:        https://github.com/fulcrumgenomics/fgsv/blob/main/LICENSE

Tools for calling breakpoints and exploring structural variation.

## Documentation

Documentation can be found in the [docs folder](docs/01_Introduction.md).

## Introduction to the `fgsv` Toolkit

The tool [`fgsv SvPileup`](https://github.com/fulcrumgenomics/fgsv/blob/main/docs/tools/SvPileup.md) takes a query-grouped BAM file as input and scans through each template one at a time, where a template is the full collection of reads and alignments from a single source molecule.
Primary and supplemental alignments for a template are used to construct a “chain” of aligned sub-segments in a way that is order and strand-aware.
These aligned sub-segments relate to each other through typical alignment mechanisms like insertions and deletions but also contain information about the relative orientation of the sub-segment to the reference genome and importantly, jumps between reference sequences (chromosomes).

For each chain of aligned sub-segments per template, outlier jumps are collected where the minimum inter-segment distance within a read must be 100bp (by default) or greater and the minimum inter-read distance per pair must be 1000bp (by default) or greater.
In the case where there is both evidence for a split-read alignment and inter-read jump, the split-read alignment evidence is favored.
At locations where these jumps occur, breakpoints are marked, and the breakpoints are given a unique ID based on the position of the breakpoint, the directionality of the left and right strands, and the other location the aligned sub-segment jumps to.
The output of this process is simply a pileup of candidate breakpoint locations.
The output of this tool is a metrics file tabulating the breakpoints and a BAM file with each alignment having custom tags that indicate which breakpoint the alignment supports (by breakpoint ID), if any.

Because of variability in short-read sequence data and their alignments, evidence for a single breakpoint may span a few loci near the true breakpoint.
The tool [`fgsv AggregateSvPileup`](https://github.com/fulcrumgenomics/fgsv/blob/main/docs/tools/AggregateSvPileup.md) is used to coalesce nearby breakpoints into one call if they appear to belong to one breakpoint.
This polishing step preserves true positive breakpoint calls and should reduce the number of false positive breakpoint calls.
Adjacent breakpoints are only merged if their left sides map to the same reference sequence, their right rides sides map to the same reference sequence, the strandedness of the left and right aligned sub-segments is the same, and their left and right positions are both within a given length threshold.
One shortcoming of the existing behavior, that should be corrected at some point, is that inter-read breakpoint evidence is considered similarly to inter-pair breakpoint evidence even though inter-read breakpoint evidence often has nucleotide-level alignment resolution and inter-pair breakpoint evidence does not.
The output of this tool is a metrics file tabulating the coalesced breakpoints with all previous breakpoint IDs listed for the new breakpoint call and an estimation of the allele frequency of the call based on the alignments that support the breakpoint.

The `fgsv` tools are an effective structural variant debugging toolkit but are not meant to be considered as a structural variant calling toolchain in-and-of-itself.
Instead, it’s better to think of the `fgsv` toolkit as an effective  “breakpoint caller”.
