package com.fulcrumgenomics.sv

import com.fulcrumgenomics.alignment.{Aligner, Alignment, Mode}
import com.fulcrumgenomics.bam.api.SamRecord
import htsjdk.samtools.CigarOperator
import htsjdk.samtools.util.SequenceUtil

/** Candidate-guided soft-clip rescue, without changing alignments or existing evidence counts. */
object JunctionRescuer {
  /** A junction alignment. Coordinates are 1-based inclusive on the selected query orientation and
    * synthetic target. Anchors count canonical query/target base columns, not deletions or insertions.
    * `reverseComplemented` indicates that the SAM-stored bases were reverse complemented.
    */
  case class Result(score: Int,
                    referenceScore: Int,
                    leftAnchor: Int,
                    rightAnchor: Int,
                    cigar: String,
                    queryStart: Int,
                    queryEnd: Int,
                    targetStart: Int,
                    targetEnd: Int,
                    reverseComplemented: Boolean,
                    identity: Double,
                    queryCoverage: Double)

  /** Cached oriented junction and the two physical reference neighborhoods. `boundary` is the actual
    * length of the incoming flank; target bases 1 through boundary belong to the left breakend.
    */
  private[sv] case class Targets(junction: Array[Byte], boundary: Int, references: IndexedSeq[Array[Byte]])

  /** Builds targets using 1-based inclusive reference intervals, clamping at contig boundaries.
    * Reference alternatives include both sides of each cut and the entire boundary-slop neighborhood.
    * The callback must return exactly the requested bases; returned arrays are not modified.
    */
  private[sv] def buildTargets(breakpoint: Breakpoint,
                   reference: (Int, Int, Int) => Array[Byte],
                   contigLengths: IndexedSeq[Int],
                   flank: Int,
                   slop: Int = 5): Targets = {
    require(breakpoint != null && reference != null && contigLengths != null, "Inputs must not be null")
    require(contigLengths.nonEmpty && contigLengths.forall(_ > 0), "Contig lengths must be positive")
    require(flank > 0 && slop >= 0, "Flank must be positive and slop non-negative")
    val ends = IndexedSeq((breakpoint.leftRefIndex, breakpoint.leftPos), (breakpoint.rightRefIndex, breakpoint.rightPos))
    ends.foreach { case (index, pos) =>
      require(index >= 0 && index < contigLengths.length, s"Invalid reference index: $index")
      require(pos >= 1 && pos <= contigLengths(index), s"Invalid position $pos on reference $index")
    }

    def fetch(index: Int, start: Long, end: Long): Array[Byte] = {
      val from = math.max(1L, start).toInt
      val to   = math.min(contigLengths(index).toLong, end).toInt
      val bases = reference(index, from, to)
      require(bases != null && bases.length == to - from + 1, s"Reference returned wrong length for $index:$from-$to")
      bases.map(b => SequenceUtil.upperCase(b))
    }

    def incoming(index: Int, pos: Int, positive: Boolean): Array[Byte] = {
      if (positive) fetch(index, pos.toLong - flank + 1, pos.toLong)
      else reverseComplement(fetch(index, pos.toLong, pos.toLong + flank - 1))
    }

    def outgoing(index: Int, pos: Int, positive: Boolean): Array[Byte] = {
      if (positive) fetch(index, pos.toLong, pos.toLong + flank - 1)
      else reverseComplement(fetch(index, pos.toLong - flank + 1, pos.toLong))
    }

    val left = incoming(breakpoint.leftRefIndex, breakpoint.leftPos, breakpoint.leftPositive)
    val right = outgoing(breakpoint.rightRefIndex, breakpoint.rightPos, breakpoint.rightPositive)
    val alternatives = ends.map { case (index, pos) =>
      fetch(index, pos.toLong - flank - slop, pos.toLong + flank + slop)
    }
    Targets(left ++ right, left.length, alternatives)
  }

  private def canonical(base: Byte): Boolean = base == 'A' || base == 'C' || base == 'G' || base == 'T'

  /** Reverse complements canonical bases in one allocation, replacing all other bases with N. */
  private def reverseComplement(bases: Array[Byte]): Array[Byte] = Array.tabulate(bases.length) { i =>
    bases(bases.length - i - 1) match {
      case 'A' => 'T'.toByte
      case 'C' => 'G'.toByte
      case 'G' => 'C'.toByte
      case 'T' => 'A'.toByte
      case _   => 'N'.toByte
    }
  }
}

/** Locally aligns whole observed reads to one candidate junction and both local reference alternatives.
  *
  * Only mapped, passing-filter, nonduplicate primary reads with an appropriate terminal soft clip are
  * considered. SAM-stored bases and their reverse complement are both tested, independently of flag 0x10.
  * Eligible read length plus slop must not exceed the configured flank, even when targets are clamped at
  * contig edges. Targets are fetched once per instance; SAM SEQ '=' additionally fetches the original
  * CIGAR-aligned reference blocks containing '=' before testing either query orientation. The reference
  * callback accepts 1-based inclusive intervals on any original alignment and must return exactly the
  * requested number of bases without restricting requests to the candidate neighborhoods. Literal reads
  * require no additional fetches. Thread safety is not provided. Local alignment
  * scores are +2 for canonical matches, -3 for mismatches, -6 for gap opening and -1 for gap extension.
  *
  * Equal best query orientations or distinct optimal alignment endpoints reject. fgbio returns only one
  * traceback per endpoint, so alternative tied gap placements at the same endpoint are not exhaustively
  * enumerated. This is evidence relative to the compared targets, not a claim of genome-wide uniqueness.
  */
class JunctionRescuer(breakpoint: Breakpoint,
                      reference: (Int, Int, Int) => Array[Byte],
                      contigLengths: IndexedSeq[Int],
                      flank: Int,
                      minAnchor: Int = 10,
                      slop: Int = 5,
                      minMapq: Int = 30,
                      minScoreDelta: Int = 10,
                      minIdentity: Double = 0.9,
                      minQueryCoverage: Double = 0.9) {
  import JunctionRescuer._

  require(minAnchor > 0 && flank >= minAnchor, "Flank must be at least the positive minimum anchor")
  require(minMapq >= 0 && minMapq <= 255, "Minimum MAPQ must be in [0,255]")
  require(minScoreDelta >= 0, "Minimum score delta must be non-negative")
  require(minIdentity > 0 && minIdentity <= 1, "Minimum identity must be in (0,1]")
  require(minQueryCoverage > 0 && minQueryCoverage <= 1, "Minimum query coverage must be in (0,1]")

  private val targets = buildTargets(breakpoint, reference, contigLengths, flank, slop)
  // Ambiguous bases must never earn a match score, including N opposite N.
  private val aligner = new Aligner(new Aligner.AlignmentScorer {
    override def scorePairing(queryBase: Byte, targetBase: Byte): Int = {
      if (canonical(queryBase) && queryBase == targetBase) 2 else -3
    }
    override def scoreGap(query: Array[Byte], target: Array[Byte], qOffset: Int, tOffset: Int,
                          gapIsInQuery: Boolean, extend: Boolean): Int = {
      if (extend) -1 else -6
    }
  }, mode=Mode.Local) {
    override protected def isMatch(queryBase: Byte, targetBase: Byte): Boolean = {
      canonical(queryBase) && queryBase == targetBase
    }
  }

  /** The actual incoming flank length after contig-edge clamping; junction target coordinates are 1-based. */
  def junctionBoundary: Int = targets.boundary

  /** Returns a rescue only when the best whole-read junction alignment passes all evidence filters.
    * Never mutates the record. Throws if an otherwise eligible read's length plus slop exceeds flank,
    * its SEQ contains '=' in an insertion/soft clip, or a reference fetch returns null/the wrong length.
    */
  def rescue(record: SamRecord): Option[Result] = {
    if (!eligible(record)) return None

    val requiredFlank = record.bases.length.toLong + slop.toLong
    require(requiredFlank <= flank.toLong,
      s"Read length ${record.bases.length} plus slop $slop exceeds flank $flank; increase --flank to at least $requiredFlank")
    val query = expandedQuery(record)
    val queries = IndexedSeq(query, reverseComplement(query))
    val best = queries.map(q => aligner.align(q, targets.junction))
    val score = best.map(_.score).max
    val winners = best.indices.filter(i => best(i).score == score)
    if (score <= 0 || winners.length != 1) return None

    val orientation = winners.head
    val alignment = best(orientation)
    val evidence = summarize(alignment)
    if (evidence.leftAnchor < minAnchor || evidence.rightAnchor < minAnchor ||
      evidence.identity < minIdentity || evidence.queryCoverage < minQueryCoverage) return None

    val referenceScores = queries.iterator.flatMap { q =>
      targets.references.iterator.map(t => aligner.align(q, t).score)
    }
    var referenceScore = 0
    while (referenceScores.hasNext) {
      referenceScore = math.max(referenceScore, referenceScores.next())
      // A tie is never evidence, even when minScoreDelta is configured as zero.
      if (score <= referenceScore || score.toLong - referenceScore.toLong < minScoreDelta.toLong) return None
    }

    // Reject separate equally optimal endpoints, even if only one of them spans the junction.
    val optima = aligner.align(queries(orientation), targets.junction, score)
    val distinct = optima.iterator.map(a => (a.queryStart, a.queryEnd, a.targetStart, a.targetEnd, a.cigar.toString)).toSet
    if (distinct.size != 1) return None

    Some(Result(score, referenceScore, evidence.leftAnchor, evidence.rightAnchor, alignment.cigar.toString,
      alignment.queryStart, alignment.queryEnd, alignment.targetStart, alignment.targetEnd, orientation == 1,
      evidence.identity, evidence.queryCoverage))
  }

  /** Resolves SAM reference-equal bases at their original aligned positions, never at the rescue target. */
  private def expandedQuery(record: SamRecord): Array[Byte] = {
    val query = record.bases.map(b => SequenceUtil.upperCase(b))
    if (query.contains('='.toByte)) {
      var queryOffset = 0
      var referencePos = record.start.toLong
      record.cigar.elems.foreach { elem =>
        val op = elem.operator
        if (op.consumesReadBases()) {
          val end = queryOffset + elem.length
          if ((queryOffset until end).exists(i => query(i) == '=')) {
            require(op.isAlignment,
              s"Read ${record.name} has SAM SEQ '=' in CIGAR operator $op; " +
                "replace '=' in soft clips or insertions with literal bases before rescue")
            val from = referencePos.toInt
            val to = (referencePos + elem.length - 1).toInt
            val bases = reference(record.refIndex, from, to)
            require(bases != null && bases.length == elem.length,
              s"Reference returned wrong length for ${record.refIndex}:$from-$to")
            (queryOffset until end).foreach { i =>
              if (query(i) == '=') query(i) = SequenceUtil.upperCase(bases(i - queryOffset))
            }
          }
          queryOffset = end
        }
        if (op.consumesReferenceBases()) referencePos += elem.length
      }
    }
    query
  }

  /** Tests observed terminal clipping at the physical cut, not the sequencing-direction strand flag. */
  private def eligible(record: SamRecord): Boolean = {
    if (!record.mapped || record.secondary || record.supplementary || record.duplicate || !record.pf ||
      record.mapq < minMapq || record.bases.isEmpty || record.start < 1 || record.refIndex < 0 ||
      record.refIndex >= contigLengths.length || record.end > contigLengths(record.refIndex)) false
    else {
      val elems = record.cigar.coalesce.elems
      val observed = elems.dropWhile(_.operator == CigarOperator.H).reverse.dropWhile(_.operator == CigarOperator.H).reverse
      val leading = observed.headOption.filter(_.operator == CigarOperator.S).map(_.length).getOrElse(0)
      val trailing = observed.lastOption.filter(_.operator == CigarOperator.S).map(_.length).getOrElse(0)
      val interior = observed.drop(if (leading > 0) 1 else 0).dropRight(if (trailing > 0) 1 else 0)
      val valid = record.cigar.elems.forall(_.length > 0) && record.cigar.lengthOnQuery == record.bases.length &&
        interior.forall(e => !e.operator.isClipping)
      // Require observed support immediately next to the clip; distant matches across a gap do not qualify.
      val leadingSource = interior.iterator.filterNot(_.operator == CigarOperator.P)
        .takeWhile(_.operator.isAlignment).map(_.length.toLong).sum
      val trailingSource = interior.reverseIterator.filterNot(_.operator == CigarOperator.P)
        .takeWhile(_.operator.isAlignment).map(_.length.toLong).sum

      def at(index: Int, pos: Int, leadingClip: Boolean): Boolean = {
        val boundary = if (leadingClip) record.start else record.end
        val clip = if (leadingClip) leading else trailing
        val source = if (leadingClip) leadingSource else trailingSource
        record.refIndex == index && math.abs(boundary.toLong - pos.toLong) <= slop &&
          clip >= minAnchor && source >= minAnchor
      }
      valid && (
        at(breakpoint.leftRefIndex, breakpoint.leftPos, leadingClip=breakpoint.leftNegative) ||
        at(breakpoint.rightRefIndex, breakpoint.rightPos, leadingClip=breakpoint.rightPositive))
    }
  }

  /** Junction evidence independent of the reference alternatives. */
  private case class Evidence(leftAnchor: Int, rightAnchor: Int, identity: Double, queryCoverage: Double)

  /** Measures aligned base columns directly so that gaps and ambiguous bases cannot inflate anchors. */
  private def summarize(alignment: Alignment): Evidence = {
    var queryOffset = alignment.queryStart - 1
    var targetOffset = alignment.targetStart - 1
    var left = 0
    var right = 0
    var matches = 0
    var columns = 0
    alignment.cigar.elems.foreach { elem =>
      val op = elem.operator
      if (op.consumesReadBases() && op.consumesReferenceBases()) {
        (0 until elem.length).foreach { offset =>
          val q = alignment.query(queryOffset + offset)
          val t = alignment.target(targetOffset + offset)
          if (canonical(q) && canonical(t)) {
            if (targetOffset + offset < targets.boundary) left += 1 else right += 1
            if (q == t) matches += 1
          }
        }
      }
      columns += elem.length
      if (op.consumesReadBases()) queryOffset += elem.length
      if (op.consumesReferenceBases()) targetOffset += elem.length
    }
    Evidence(left, right,
      if (columns == 0) 0.0 else matches.toDouble / columns,
      alignment.cigar.lengthOnQuery.toDouble / alignment.query.length)
  }
}
