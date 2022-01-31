package com.fulcrumgenomics.fgsv.tools

import com.fulcrumgenomics.FgBioDef.{PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.fgsv.cmdline.{ClpGroups, SVTool}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}

@clp(group=ClpGroups.All, description=
  """
    |Moves the UMI at the end of the read name to the RX tag.
  """)
class MoveUmiFromReadName
( @arg(flag='i', doc="The input BAM file") input: PathToBam,
  @arg(flag='o', doc="The output BAM file") output: PathToBam,
) extends SVTool {
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val source = SamSource(input)
    val writer = SamWriter(output, source.header)
    val progress = new ProgressLogger(logger)
    source.foreach { rec =>
      progress.record(rec)
      val fields = rec.name.split(':')
      require(fields.length > 1)
      rec.name  = fields.iterator.take(fields.length - 1).mkString(":")
      rec("RX") = fields.last
      writer += rec
    }
    source.safelyClose()
    writer.close()
  }
}