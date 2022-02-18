# Contributing

<!--toc start-->
- [Introduction](#introduction)
- [Updating the Metrics Markdown](#updating-the-metrics-markdown)
- [Updating the Tools Markdown](#updating-the-tools-markdown)
<!--toc end-->

## Introduction

This code for this repository is written in [Scala][scala-link] leveraging the [fgbio][fgbio-link] toolkit.

All contributions must be tested, documented, and code reviewed. 

No changes should be pushed to the `main` branch.

When the definition or documentation for metrics or command line tools are updated or added, please update the 
Markdown documentation below for [metrics](#updating-the-metrics-markdown) and [tools](#updating-the-tools-markdown).

## Updating the Metrics Markdown

The [metric descriptions listed here][metrics-link] is auto-generated with the following script:

```
bash src/scripts/build_metric_docs.sh
```

You will need to install `jq` via `sudo apt-get install jq` (for Debian/Ubuntu) or `brew install jq` for OSX.

## Updating the Tools Markdown

The [tool descriptions listed here][tools-link] is auto-generated with the following script:

```
bash src/scripts/build_tool_docs.sh
```

[metrics-link]: 04_Metrics.md
[tools-link]:   05_Tools.md
[scala-link]:   https://www.scala-lang.org/
[fgbio-link]:   https://github.com/fulcrumgenomics/fgbio/
