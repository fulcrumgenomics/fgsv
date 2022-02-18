/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package com.fulcrumgenomics.sv.internal

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.Sopt
import com.fulcrumgenomics.sopt.Sopt.{CommandSuccess, Failure}

/** Trait for internal tools to extends. */
trait InternalTool extends LazyLogging {
  def execute(): Unit
}

/** Main class for internal tools. */
object FgSvInternalMain {
  def main(args: Array[String]): Unit = {
    val tools = Sopt.find[InternalTool](packages=Seq("com.fulcrumgenomics.sv.internal"))
    Sopt.parseCommand[InternalTool]("fgsv-internal", args.toIndexedSeq, tools) match {
      case Failure(usage) => {
        System.err.println(usage())
        System.exit(1)
      }
      case CommandSuccess(tool) =>
        tool.execute()
      case _ => unreachable()
    }
  }
}