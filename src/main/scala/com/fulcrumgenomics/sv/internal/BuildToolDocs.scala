package com.fulcrumgenomics.sv.internal

import com.fulcrumgenomics.FgBioDef.DirPath
import com.fulcrumgenomics.internal.{BuildToolDocs => FgBioBuildToolDocs}
import com.fulcrumgenomics.sopt.{arg, clp}

@clp(description="Generates the suite of per-tool MarkDown documents.")
class BuildToolDocs
( @arg(flag='o', doc="Output directory") output: DirPath,
  @arg(flag='p', doc="The packages to document") packages: Seq[String] = Seq("com.fulcrumgenomics.sv"),
  @arg(flag='n', doc="The name of the tool chain") name: String = "fgsv"
) extends FgBioBuildToolDocs(output=output, packages=packages, name=name) with InternalTool
