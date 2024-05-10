# fgsv

[![Bioconda][bioconda-badge-link]][bioconda-link]
[![Build Status][github-badge]][github-link]
[![Language][scala-badge]][scala-link]
[![License][license-badge]][license-link]
[![DOI][doi-badge]][doi-link]

[bioconda-badge-link]: https://img.shields.io/conda/dn/bioconda/fgsv.svg?label=Bioconda
[bioconda-link]:       http://bioconda.github.io/recipes/fgsv/README.html
[github-badge]:        https://github.com/fulcrumgenomics/fgsv/actions/workflows/unittests.yaml/badge.svg?branch=main
[github-link]:         https://github.com/fulcrumgenomics/fgsv/actions/workflows/unittests.yaml
[scala-badge]:         https://img.shields.io/badge/language-scala-c22d40.svg
[scala-link]:          https://www.scala-lang.org/
[license-badge]:       https://img.shields.io/badge/license-MIT-blue.svg
[license-link]:        https://github.com/fulcrumgenomics/fgsv/blob/main/LICENSE
[doi-badge]:           https://zenodo.org/badge/454071954.svg
[doi-link]:            https://zenodo.org/doi/10.5281/zenodo.10452647

Tools to gather evidence for structural variation via breakpoint detection.

## Documentation

Documentation can be found in the [docs folder](docs/01_Introduction.md).

## Introduction to the `fgsv` Toolkit

The `fgsv` toolkit contains tools for effective structural variant debugging but are not meant to be used as a structural variant calling toolchain in-and-of-itself.
Instead, it is better to think of `fgsv` as an effective breakpoint detection and structural variant exploration toolkit.

When describing structural variation, we use the term breakpoint to mean a junction between two loci and the term breakend to refer to one of the loci in a breakpoint.
Importantly, all point intervals (1-length) reported by this toolkit are 1-based inclusive from the perspective of the reference sequence.

### `fgsv SvPileup`

Collates a pileup of putative structural variant supporting reads.

```console
fgsv SvPileup \
    --input sample.bam \
    --output sample.svpileup
```

The tool [`fgsv SvPileup`](https://github.com/fulcrumgenomics/fgsv/blob/main/docs/tools/SvPileup.md) takes a query-grouped BAM file as input and scans through each template one at a time, where a template is the full collection of reads and alignments from a single source molecule.
For example, a paired-end read may have an alignment per read: one alignment for read 1 and another alignment for read 2.

Primary and supplementary alignments for a template (see the [SAM Format Specification v1](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information) are used to construct a “chain” of aligned sub-segments in a way that honors the logical ordering of sub-segments and their strandeness in relation to the reference sequence.
These aligned sub-segments in a chain relate to each other through typical alignment mechanisms like insertions and deletions but also contain information about the relative orientation of the sub-segment to the reference sequence and importantly, jumps between reference sequences such as translocations between chromosomes or contigs.

For each chain of aligned sub-segments per template, outlier jumps are collected where the minimum inter-segment distance within a read must be 100bp (by default) or greater, and the minimum inter-read distance across reads (e.g. between reads in a paired-end read) must be 1000bp (by default) or greater.
In the case where there is both evidence for a split-read alignment and inter-read jump, the split-read alignment evidence is favored since it gives a precise breakpoint.
At locations where these jumps occur, breakpoints are marked and the breakpoints are given a unique ID based on the positions of the breakends and the directionality of the left and right strands leading into each breakend.

This process creates a collection of candidate breakpoint locations.
The output of this tool is a metrics file tabulating the breakpoints and a BAM file with each breakpoint-supporting alignment having custom tags that indicate which breakpoint the alignment supports.

### `fgsv AggregateSvPileup`

Merges nearby pileups of reads supporting putative breakpoints.

```console
fgsv AggregateSvPileup \
    --bam sample.bam \
    --input sample.svpileup.txt \
    --output sample.svpileup.aggregate.txt
```

Because of variability in typical short-read alignments, evidence for a single breakpoint may span a few loci near the true breakend loci. For example, if the breakpoint only has intra-read evidence, then the breakpoint could coincidentally occur within the unobserved bases between read 1 and read 2 in a pair. In other cases and due to sequence similarity or homology between each breakend locus, it is not always possible to locate the exact nucleotide point where the breakends occur, and instead a plausible region may exist that supports either breakend loci.

The tool [`fgsv AggregateSvPileup`](https://github.com/fulcrumgenomics/fgsv/blob/main/docs/tools/AggregateSvPileup.md) is used to coalesce nearby breakpoints into one event if they appear to belong to one true breakpoint.
This polishing step preserves true positive breakpoint events and intends to reduce the number of false positive breakpoint events.

Adjacent breakpoints are only merged if their left breakends map to the same reference sequence, their right breakends map to the same reference sequence, the strandedness of the left and right aligned sub-segments is the same, and their left and right genomic breakend positions are both within a given length threshold.

One shortcoming of the existing behavior, which should be corrected at some point, is that intra-read breakpoint evidence is considered similarly to inter-pair breakpoint evidence even though intra-read breakpoint evidence often has nucleotide-level alignment resolution and inter-pair breakpoint evidence does not.

The output of this tool is a metrics file tabulating the coalesced breakpoints with all previous breakpoint IDs listed for the new breakpoint event and an estimation of the allele frequency of the event based on the alignments that support the breakpoint.

## `AggregateSvPileupToBedPE`

Convert the output of `AggregateSvPileup` to BEDPE.

```console
fgsv AggregateSvPileupToBedPE \
    --input sample.svpileup.aggregate.txt \
    --output sample.svpileup.aggregate.bedpe
```

The tool [`fgsv AggregateSvPileupToBedPE`](https://github.com/fulcrumgenomics/fgsv/blob/main/docs/tools/AggregateSvPileupToBedPE.md) is used to convert the output of `AggregateSvPileup` to the [BEDPE format](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) so that it can be viewed in [IGV](https://igv.org/) and other BEDPE-supporting genome browsers. For example:

![BEDPE in IGV](docs/img/fgsv-bedpe.png)
