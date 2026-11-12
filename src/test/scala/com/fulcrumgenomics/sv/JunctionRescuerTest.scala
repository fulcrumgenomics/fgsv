package com.fulcrumgenomics.sv

import com.fulcrumgenomics.alignment.{Aligner, Cigar, Mode}
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.testing.SamBuilder
import com.fulcrumgenomics.util.Sequences

import scala.util.Random

class JunctionRescuerTest extends UnitSpec {
  private def reverseComplement(bases: String): String = Sequences.revcomp(bases)

  // Independent pinned fgbio aligner for canonical fixture preconditions, not the rescue implementation.
  private def fixtureAligner: Aligner = Aligner(2, -3, -6, -1, Mode.Local)

  private def randomBases(length: Int, seed: Int): String = {
    val random = new Random(seed.toLong)
    (0 until length).map(_ => "ACGT"(random.nextInt(4))).mkString
  }

  /** Generates sequence independently of buildTargets, and records in SAM reference orientation. */
  private class Fixture(val bp: Breakpoint = Breakpoint(0, 600, true, 1, 1400, true), seed: Int = 17) {
    val refs: Array[String] = Array(randomBases(2500, seed), randomBases(2500, seed + 1))
    val lengths: IndexedSeq[Int] = refs.map(_.length).toIndexedSeq
    val reference: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
      require(start >= 1 && end <= refs(index).length && end >= start)
      refs(index).substring(start - 1, end).getBytes("US-ASCII")
    }

    def rescuer: JunctionRescuer = new JunctionRescuer(bp, reference, lengths, flank=200, minAnchor=20)

    def left(length: Int): String = {
      if (bp.leftPositive) refs(bp.leftRefIndex).substring(bp.leftPos - length, bp.leftPos)
      else reverseComplement(refs(bp.leftRefIndex).substring(bp.leftPos - 1, bp.leftPos - 1 + length))
    }

    def right(length: Int): String = {
      if (bp.rightPositive) refs(bp.rightRefIndex).substring(bp.rightPos - 1, bp.rightPos - 1 + length)
      else reverseComplement(refs(bp.rightRefIndex).substring(bp.rightPos - length, bp.rightPos))
    }

    /** A contiguous reference read whose terminal bases were soft clipped at the supplied cut. */
    def referenceRead(sourceLeft: Boolean, clipLength: Int = 10, offset: Int = 0): SamRecord = {
      val record = sourceRead(sourceLeft, clipLength=clipLength)
      record.start += offset
      val leadingClip = record.cigar.elems.headOption.filter(_.operator == htsjdk.samtools.CigarOperator.S)
        .map(_.length).getOrElse(0)
      val from = record.start - leadingClip - 1
      record.bases = refs(record.refIndex).substring(from, from + record.bases.length)
      record
    }

    /** Places the source anchor and opposite clip on either candidate side. */
    def sourceRead(sourceLeft: Boolean, sourceLength: Int = 125, clipLength: Int = 25,
                   negativeFlag: Boolean = false): SamRecord = {
      read(leftLength=if (sourceLeft) sourceLength else clipLength,
        rightLength=if (sourceLeft) clipLength else sourceLength, sourceLeft=sourceLeft, negativeFlag=negativeFlag)
    }

    def read(leftLength: Int = 125, rightLength: Int = 25, sourceLeft: Boolean = true,
             negativeFlag: Boolean = false): SamRecord = {
      val positive = if (sourceLeft) bp.leftPositive else bp.rightPositive
      val sequence = left(leftLength) + right(rightLength)
      val stored = if (positive) sequence else reverseComplement(sequence)
      val source = if (sourceLeft) leftLength else rightLength
      val clip = if (sourceLeft) rightLength else leftLength
      val leading = if (sourceLeft) !positive else positive
      val pos = if (sourceLeft) bp.leftPos else bp.rightPos
      val cigar = if (leading) s"${clip}S${source}M" else s"${source}M${clip}S"
      val builder = new SamBuilder(readLength=stored.length)
      builder.addFrag(bases=stored, quals="I" * stored.length,
        contig=if (sourceLeft) bp.leftRefIndex else bp.rightRefIndex,
        start=if (leading) pos else pos - source + 1, cigar=cigar,
        strand=if (negativeFlag) SamBuilder.Minus else SamBuilder.Plus).value
    }
  }

  "JunctionRescuer" should "rescue all strand pairs, both candidate sides, and both SAM strand flags" in {
    for (leftPositive <- IndexedSeq(true, false);
         rightPositive <- IndexedSeq(true, false);
         sourceLeft <- IndexedSeq(true, false);
         negativeFlag <- IndexedSeq(true, false)) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive))
      val record = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=25, negativeFlag=negativeFlag)
      val before = record.asSam.getSAMString
      val result = fixture.rescuer.rescue(record).value
      withClue(s"$leftPositive/$rightPositive sourceLeft=$sourceLeft flag=$negativeFlag: ") {
        result.score shouldBe 300
        result.leftAnchor shouldBe (if (sourceLeft) 125 else 25)
        result.rightAnchor shouldBe (if (sourceLeft) 25 else 125)
        result.reverseComplemented shouldBe !(if (sourceLeft) leftPositive else rightPositive)
        result.queryStart shouldBe 1
        result.queryEnd shouldBe 150
        result.targetStart shouldBe 201 - result.leftAnchor
        result.targetEnd shouldBe 200 + result.rightAnchor
        result.identity shouldBe 1.0
        result.queryCoverage shouldBe 1.0
        result.referenceScore should be <= 290
        record.asSam.getSAMString shouldBe before
      }
    }
  }

  it should "resolve reference-equal SEQ with padded and split-clip CIGARs without changing inputs" in {
    for (leftPositive <- IndexedSeq(true, false); rightPositive <- IndexedSeq(true, false);
         sourceLeft <- IndexedSeq(true, false); negativeFlag <- IndexedSeq(true, false)) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive))
      val original = fixture.sourceRead(sourceLeft, negativeFlag=negativeFlag)
      val expected = fixture.rescuer.rescue(original).value
      val leading = original.cigar.toString.startsWith("25S")
      for (split <- IndexedSeq(true, false); padded <- IndexedSeq(true, false)) {
        val clip = if (split) "12S13S" else "25S"
        val source = if (!padded) "40M45=40X" else if (leading) "2P10M2P35=80X" else "80X35=2P10M2P"
        val cigar = "3H" + (if (leading) clip + source else source + clip) + "4H"
        for (equal <- IndexedSeq(true, false)) {
          val record = fixture.sourceRead(sourceLeft, negativeFlag=negativeFlag)
          record.cigar = Cigar(cigar)
          if (equal) record.bases = if (leading) record.basesString.take(25) + ("=" * 125)
            else ("=" * 125) + record.basesString.takeRight(25)
          val before = record.asSam.getSAMString
          val fetched = scala.collection.mutable.ArrayBuffer.empty[(Array[Byte], Array[Byte])]
          var calls = 0
          val callback: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
            calls += 1
            val bases = fixture.reference(index, start, end).map(b => (b + 32).toByte)
            fetched += ((bases, bases.clone()))
            bases
          }
          val rescuer = new JunctionRescuer(fixture.bp, callback, fixture.lengths, flank=200, minAnchor=20)
          withClue(s"strands=$leftPositive/$rightPositive sourceLeft=$sourceLeft flag=$negativeFlag cigar=$cigar equal=$equal: ") {
            rescuer.rescue(record).value shouldBe expected
            calls shouldBe (if (equal) 7 else 4)
            record.asSam.getSAMString shouldBe before
            fetched.foreach { case (bases, snapshot) => bases.toIndexedSeq shouldBe snapshot.toIndexedSeq }
          }
        }
      }
    }
  }

  it should "expand reference-equal bases at original coordinates across insertions deletions skips and padding" in {
    for (sourceLeft <- IndexedSeq(true, false); negativeFlag <- IndexedSeq(true, false);
         gap <- IndexedSeq("2I", "2D", "2N", "300N")) {
      val fixture = new Fixture()
      val record = fixture.sourceRead(sourceLeft, negativeFlag=negativeFlag)
      record.cigar = Cigar(if (sourceLeft) s"100M${gap}2P25M25S" else s"25S25M2P${gap}100M")
      if (sourceLeft) record.start = fixture.bp.leftPos - record.cigar.lengthOnTarget + 1
      val intervals = scala.collection.mutable.ArrayBuffer.empty[(Int, Int, Int)]
      val alignedIntervals = scala.collection.mutable.ArrayBuffer.empty[(Int, Int, Int)]
      val literal = new StringBuilder
      val encoded = new StringBuilder
      var pos = record.start
      record.cigar.elems.foreach { elem =>
        val op = elem.operator
        if (op.isAlignment) {
          alignedIntervals += ((record.refIndex, pos, pos + elem.length - 1))
          literal.append(new String(fixture.reference(record.refIndex, pos, pos + elem.length - 1), "US-ASCII"))
          encoded.append("=" * elem.length)
        }
        else if (op == htsjdk.samtools.CigarOperator.S) {
          val clip = if (sourceLeft) record.basesString.takeRight(25) else record.basesString.take(25)
          literal.append(clip)
          encoded.append(clip)
        }
        else if (op.consumesReadBases()) {
          literal.append("AC")
          encoded.append("AC")
        }
        if (op.consumesReferenceBases()) pos += elem.length
      }
      record.bases = literal.toString
      record.quals = "I" * literal.length
      val expected = fixture.rescuer.rescue(record)
      if (gap != "300N") expected.isDefined shouldBe true
      val callback: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
        intervals += ((index, start, end))
        fixture.reference(index, start, end)
      }
      val rescuer = new JunctionRescuer(fixture.bp, callback, fixture.lengths, flank=200, minAnchor=20)
      intervals.clear()
      record.bases = encoded.toString
      val before = record.asSam.getSAMString
      withClue(s"sourceLeft=$sourceLeft flag=$negativeFlag gap=$gap: ") {
        rescuer.rescue(record) shouldBe expected
        intervals.toIndexedSeq shouldBe alignedIntervals.toIndexedSeq
        record.asSam.getSAMString shouldBe before
      }
    }
  }

  it should "preserve literal mismatches and fetch only aligned blocks containing reference-equal bases" in {
    val fixture = new Fixture()
    for (sourceLeft <- IndexedSeq(true, false)) {
      val record = fixture.sourceRead(sourceLeft)
      record.cigar = Cigar(if (sourceLeft) "50M25=50X25S" else "25S50M25=50X")
      val sourceOffset = if (sourceLeft) 0 else 25
      val mismatchOffset = sourceOffset + 10
      val literal = record.basesString
      record.bases = literal.updated(mismatchOffset, if (literal(mismatchOffset) == 'A') 'C' else 'A')
      val expected = fixture.rescuer.rescue(record).value
      expected.score shouldBe 295
      val encoded = record.basesString.zipWithIndex.map { case (base, i) =>
        if (i >= sourceOffset + 50 && i < sourceOffset + 75 && i % 2 == 0) '=' else base
      }.mkString
      val intervals = scala.collection.mutable.ArrayBuffer.empty[(Int, Int, Int)]
      val callback: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
        intervals += ((index, start, end))
        fixture.reference(index, start, end)
      }
      val rescuer = new JunctionRescuer(fixture.bp, callback, fixture.lengths, flank=200, minAnchor=20)
      intervals.clear()
      record.bases = encoded
      val before = record.asSam.getSAMString
      rescuer.rescue(record).value shouldBe expected
      intervals.toIndexedSeq shouldBe IndexedSeq((record.refIndex, record.start + 50, record.start + 74))
      record.asSam.getSAMString shouldBe before
    }
  }

  it should "reject reference-equal SEQ in clips and insertions with actionable errors" in {
    val fixture = new Fixture()
    for (sourceLeft <- IndexedSeq(true, false); insertion <- IndexedSeq(true, false)) {
      val record = fixture.sourceRead(sourceLeft)
      val offset = if (insertion) {
        record.cigar = Cigar(if (sourceLeft) "99M1I25M25S" else "25S25M1I99M")
        if (sourceLeft) record.start = fixture.bp.leftPos - record.cigar.lengthOnTarget + 1
        if (sourceLeft) 99 else 50
      }
      else if (sourceLeft) 149 else 0
      record.bases = record.basesString.updated(offset, '=')
      val before = record.asSam.getSAMString
      val error = intercept[IllegalArgumentException](fixture.rescuer.rescue(record))
      error.getMessage should include (s"Read ${record.name}")
      error.getMessage should include (if (insertion) "operator I" else "operator S")
      error.getMessage should include ("literal bases before rescue")
      record.asSam.getSAMString shouldBe before
    }
  }

  it should "validate reference lengths when expanding reference-equal SEQ" in {
    val fixture = new Fixture()
    for (returned <- IndexedSeq(null, Array.emptyByteArray, Array.fill[Byte](124)('A'.toByte))) {
      val record = fixture.read()
      record.bases = ("=" * 125) + record.basesString.takeRight(25)
      var calls = 0
      val callback: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
        calls += 1
        if (calls <= 4) fixture.reference(index, start, end) else returned
      }
      val rescuer = new JunctionRescuer(fixture.bp, callback, fixture.lengths, flank=200)
      val before = record.asSam.getSAMString
      val error = intercept[IllegalArgumentException](rescuer.rescue(record))
      error.getMessage should include ("Reference returned wrong length for 0:476-600")
      record.asSam.getSAMString shouldBe before
    }
  }

  it should "rescue generated DEL, DUP, TRA and inversion reads in either canonical representation" in {
    val breakpoints = IndexedSeq(
      Breakpoint(0, 600, true, 0, 1400, true),
      Breakpoint(0, 1400, true, 0, 600, true),
      Breakpoint(0, 600, true, 1, 1400, true),
      Breakpoint(0, 600, true, 0, 1400, false),
      Breakpoint(0, 600, false, 0, 1400, true),
      Breakpoint(0, 600, false, 1, 1400, false))
    for (bp <- breakpoints; seed <- 1 to 4; sourceLeft <- IndexedSeq(true, false)) {
      val fixture = new Fixture(bp, seed)
      val record = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=25)
      fixture.rescuer.rescue(record).value.score shouldBe 300
      val reversed = new JunctionRescuer(bp.reversed, fixture.reference, fixture.lengths, flank=200)
      val result = reversed.rescue(record).value
      result.score shouldBe 300
      result.leftAnchor shouldBe (if (sourceLeft) 25 else 125)
      result.rightAnchor shouldBe (if (sourceLeft) 125 else 25)
    }
  }

  it should "default to ten aligned bases on each side" in {
    val fixture = new Fixture()
    val rescuer = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200)
    for (sourceLeft <- IndexedSeq(true, false)) {
      val accepted = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=10)
      rescuer.rescue(accepted).isDefined shouldBe true
      val rejected = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=9)
      rescuer.rescue(rejected) shouldBe None
    }
  }

  it should "rescue the actual 25S125M geometry and exact minimum 20-base anchors" in {
    val fixture = new Fixture()
    val record = fixture.read(leftLength=25, rightLength=125, sourceLeft=false)
    record.cigar.toString shouldBe "25S125M"
    fixture.rescuer.rescue(record).value.leftAnchor shouldBe 25
    for (sourceLeft <- IndexedSeq(true, false)) {
      val read = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=20)
      fixture.rescuer.rescue(read).value.score shouldBe 290
    }
  }

  it should "reject 1, 3 and 19-base observed clips and too-short source anchors on both sides" in {
    val fixture = new Fixture()
    for (length <- IndexedSeq(1, 3, 19); sourceLeft <- IndexedSeq(true, false)) {
      val shortClip = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=length)
      fixture.rescuer.rescue(shortClip) shouldBe None
      val shortSource = fixture.sourceRead(sourceLeft, sourceLength=length, clipLength=125)
      fixture.rescuer.rescue(shortSource) shouldBe None
    }
  }

  it should "allow outer hard clips but never count them as observed sequence" in {
    val fixture = new Fixture()
    for (sourceLeft <- IndexedSeq(true, false)) {
      val read = fixture.sourceRead(sourceLeft)
      read.cigar = Cigar("12H" + read.cigar.toString + "8H")
      fixture.rescuer.rescue(read).value.score shouldBe 300
      val short = fixture.read(leftLength=125, rightLength=3)
      short.cigar = Cigar("125M3S30H")
      fixture.rescuer.rescue(short) shouldBe None
      short.cigar = Cigar("125M33H")
      short.bases = short.bases.take(125)
      fixture.rescuer.rescue(short) shouldBe None
    }
  }

  it should "reject reads that fail primary alignment, mapping, duplicate or QC gates" in {
    val fixture = new Fixture()
    val changes: IndexedSeq[SamRecord => Unit] = IndexedSeq(
      r => r.mapq = 29,
      r => r.duplicate = true,
      r => r.pf = false,
      r => r.supplementary = true,
      r => r.secondary = true,
      r => r.unmapped = true,
      r => r.refIndex = 2,
      r => r.cigar = Cigar("150M"))
    changes.foreach { change =>
      val record = fixture.read()
      change(record)
      fixture.rescuer.rescue(record) shouldBe None
    }
    val threshold = fixture.read()
    threshold.mapq = 30
    fixture.rescuer.rescue(threshold).isDefined shouldBe true
  }

  it should "gate the correct clip boundary at exactly slop and reject one base outside" in {
    for (leftPositive <- IndexedSeq(true, false); rightPositive <- IndexedSeq(true, false);
         sourceLeft <- IndexedSeq(true, false); offset <- IndexedSeq(-6, -5, 5, 6)) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive))
      val record = fixture.sourceRead(sourceLeft)
      record.start += offset
      fixture.rescuer.rescue(record).isDefined shouldBe (math.abs(offset) <= 5)
      val wrongClip = fixture.sourceRead(sourceLeft)
      wrongClip.cigar = Cigar(wrongClip.cigar.elems.reverse)
      fixture.rescuer.rescue(wrongClip) shouldBe None
    }
  }

  it should "reject wrong opposite sequence and a best short one-sided local alignment" in {
    val fixture = new Fixture()
    val record = fixture.read()
    record.bases = fixture.left(125) + ("N" * 25)
    fixture.rescuer.rescue(record) shouldBe None
    record.bases = fixture.left(125) + randomBases(25, 999)
    fixture.rescuer.rescue(record) shouldBe None
  }

  it should "reject eligible reference-only clips near the cut while rescuing ten-base positives" in {
    for (seed <- 31 to 33; leftPositive <- IndexedSeq(true, false); rightPositive <- IndexedSeq(true, false);
         sourceLeft <- IndexedSeq(true, false)) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive), seed)
      val rescuer = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200)
      for (offset <- IndexedSeq(-5, 0, 5)) {
        withClue(s"seed=$seed strands=$leftPositive/$rightPositive sourceLeft=$sourceLeft offset=$offset: ") {
          rescuer.rescue(fixture.referenceRead(sourceLeft, offset=offset)) shouldBe None
          val positive = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=10)
          positive.start += offset
          val result = rescuer.rescue(positive).value
          math.min(result.leftAnchor, result.rightAnchor) shouldBe 10
          result.score shouldBe 270
        }
      }
    }
  }

  it should "reject repeated opposite-flank sequence with an equally good local reference explanation" in {
    for (seed <- 41 to 43; clipLength <- IndexedSeq(10, 25)) {
      val fixture = new Fixture(seed=seed)
      val repeat = randomBases(5, seed + 100) * 20
      fixture.refs(0) = fixture.refs(0).patch(fixture.bp.leftPos, repeat, repeat.length)
      fixture.refs(1) = fixture.refs(1).patch(fixture.bp.rightPos - 1, repeat, repeat.length)
      val record = fixture.read(rightLength=clipLength)
      val rescuer = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minScoreDelta=0)
      rescuer.rescue(record) shouldBe None
    }
  }

  it should "prefer the whole-read reference explanation for nearby insertions and deletions" in {
    for (seed <- 51 to 53; change <- IndexedSeq("insertion", "deletion")) {
      val fixture = new Fixture(seed=seed)
      val record = fixture.referenceRead(sourceLeft=true, clipLength=25)
      val original = record.basesString
      if (change == "insertion") {
        record.bases = original.patch(100, "AC", 0)
        record.cigar = Cigar("100M2I25M25S")
      }
      else {
        record.bases = original.patch(100, "", 2)
        record.cigar = Cigar("100M2D23M25S")
      }
      new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200)
        .rescue(record) shouldBe None
    }
  }

  it should "reject adapter-like tails and erroneous or ambiguous short clips" in {
    for (seed <- 61 to 63; sourceLeft <- IndexedSeq(true, false)) {
      val fixture = new Fixture(seed=seed)
      val rescuer = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200)
      val opposite = if (sourceLeft) fixture.right(10) else fixture.left(10)
      val substitutions = opposite.zipWithIndex.map { case (base, index) =>
        if (index % 3 == 1) { if (base == 'A') 'C' else 'A' } else base
      }.mkString
      val tails = IndexedSeq("AGATCGGAAG", "N" * 10, opposite.updated(5, 'N'), substitutions)
      tails.foreach { tail =>
        val record = fixture.sourceRead(sourceLeft, sourceLength=125, clipLength=10)
        record.bases = if (sourceLeft) fixture.left(125) + tail else tail + fixture.right(125)
        withClue(s"seed=$seed sourceLeft=$sourceLeft tail=$tail: ") {
          rescuer.rescue(record) shouldBe None
        }
      }
    }
  }

  it should "require score advantage even for a perfect ten-base opposite anchor" in {
    val fixture = new Fixture(seed=71)
    val record = fixture.read(rightLength=10)
    val query = record.basesString
    val decoy = query.updated(60, if (query(60) == 'A') 'C' else 'A')
    fixture.refs(0) = fixture.refs(0).patch(fixture.bp.leftPos + 30, decoy, decoy.length)
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200)
      .rescue(record) shouldBe None
    val result = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minScoreDelta=5)
      .rescue(record).value
    result.rightAnchor shouldBe 10
    result.score shouldBe 270
    result.referenceScore shouldBe 265
  }

  it should "reject a contiguous reference tie even when the score delta is zero" in {
    val fixture = new Fixture(Breakpoint(0, 600, true, 0, 601, true))
    val record = fixture.read()
    fixture.rescuer.rescue(record) shouldBe None
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minScoreDelta=0)
      .rescue(record) shouldBe None
  }

  it should "compare both reference alternatives in both orientations" in {
    for (leftPositive <- IndexedSeq(true, false); rightPositive <- IndexedSeq(true, false);
         decoyLeft <- IndexedSeq(true, false); reversed <- IndexedSeq(true, false)) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive))
      val record = fixture.read()
      val query = record.basesString
      val decoy = if (reversed) reverseComplement(query) else query
      val junction = fixture.left(200) + fixture.right(200)
      fixture.rescuer.rescue(record).value.score shouldBe 300
      val index = if (decoyLeft) fixture.bp.leftRefIndex else fixture.bp.rightRefIndex
      val pos = if (decoyLeft) fixture.bp.leftPos else fixture.bp.rightPos
      // Put the entire decoy on the physical REF-only side, outside the ALT flank.
      val afterCut = if (decoyLeft) leftPositive else !rightPositive
      val from = if (afterCut) pos + 35 else pos - 35 - decoy.length + 1
      fixture.refs(index) = fixture.refs(index).patch(from - 1, decoy, decoy.length)
      withClue(s"strands=$leftPositive/$rightPositive decoyLeft=$decoyLeft reversed=$reversed: ") {
        fixture.left(200) + fixture.right(200) shouldBe junction
        val aligner = fixtureAligner
        val queries = IndexedSeq(query, reverseComplement(query)).map(_.getBytes("US-ASCII"))
        val altScores = queries.map(q => aligner.align(q, junction.getBytes("US-ASCII")).score)
        altScores.max shouldBe 300
        altScores.count(_ == 300) shouldBe 1
        val winner = queries(altScores.indexOf(300))
        aligner.align(winner, junction.getBytes("US-ASCII"), 300).size shouldBe 1
        val references = IndexedSeq((0, 600), (1, 1400)).map { case (refIndex, cut) =>
          fixture.reference(refIndex, cut - 205, cut + 205)
        }
        for (orientation <- queries.indices; refIndex <- references.indices) {
          val score = aligner.align(queries(orientation), references(refIndex)).score
          if (refIndex == index && orientation == (if (reversed) 1 else 0)) score shouldBe 300
          else score should be <= 290
        }
        fixture.rescuer.rescue(record) shouldBe None
      }
    }
  }

  it should "require the configured score advantage" in {
    val fixture = new Fixture()
    val record = fixture.read()
    val result = fixture.rescuer.rescue(record).value
    val delta = result.score - result.referenceScore
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minScoreDelta=delta)
      .rescue(record).isDefined shouldBe true
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minScoreDelta=delta + 1)
      .rescue(record) shouldBe None
  }

  it should "use Smith-Waterman for mismatches and internal insertions and deletions" in {
    for (change <- IndexedSeq("mismatch", "insertion", "deletion")) {
      val fixture = new Fixture()
      val record = fixture.read(leftLength=75, rightLength=75)
      val original = record.basesString
      val sequence = change match {
        case "mismatch"  => original.updated(35, if (original(35) == 'A') 'C' else 'A')
        case "insertion" => original.patch(35, "AC", 0)
        case _           => original.patch(35, "", 2)
      }
      record.bases = sequence
      record.cigar = Cigar(s"${sequence.length - 75}M75S")
      record.start = fixture.bp.leftPos - (sequence.length - 75) + 1
      val result = fixture.rescuer.rescue(record).value
      result.score should be < 300
      result.identity should be >= 0.9
      result.queryCoverage shouldBe 1.0
      if (change == "insertion") result.cigar should include ("I")
      if (change == "deletion") result.cigar should include ("D")
      if (change == "mismatch") result.identity should be < 1.0
    }
  }

  it should "count aligned base columns rather than deletion length as anchors" in {
    val fixture = new Fixture()
    val sequence = fixture.left(25).take(10) + fixture.left(25).takeRight(5) + fixture.right(125)
    val record = fixture.read(leftLength=25, rightLength=125, sourceLeft=false)
    record.bases = sequence
    // Twenty observed clipped bases, but at most fifteen can align to the incoming flank.
    record.cigar = Cigar("20S120M")
    val permissive = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200,
      minAnchor=10, minIdentity=0.5, minQueryCoverage=0.5, minScoreDelta=0).rescue(record).value
    permissive.cigar should include ("D")
    permissive.leftAnchor should be < 20
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200,
      minAnchor=20, minIdentity=0.5, minQueryCoverage=0.5, minScoreDelta=0).rescue(record) shouldBe None
    val deletedSource = fixture.read()
    deletedSource.cigar = Cigar("19M106D131S")
    fixture.rescuer.rescue(deletedSource) shouldBe None
  }

  it should "require the source anchor immediately adjacent to the soft clip" in {
    for (leftPositive <- IndexedSeq(true, false); rightPositive <- IndexedSeq(true, false);
         sourceLeft <- IndexedSeq(true, false); negativeFlag <- IndexedSeq(true, false);
         gap <- IndexedSeq("I", "D", "N")) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive))
      val sourceLength = if (gap == "I") 120 else 119
      val record = fixture.sourceRead(sourceLeft, sourceLength=sourceLength, clipLength=30, negativeFlag=negativeFlag)
      val leading = record.cigar.toString.startsWith("30S")
      record.cigar = Cigar(if (leading) s"15S15S19M2P1${gap}2P100M" else s"100M2P1${gap}2P19M15S15S")
      if (!leading) {
        val pos = if (sourceLeft) fixture.bp.leftPos else fixture.bp.rightPos
        record.start = pos - record.cigar.lengthOnTarget + 1
      }
      val before = record.asSam.getSAMString
      withClue(s"strands=$leftPositive/$rightPositive sourceLeft=$sourceLeft flag=$negativeFlag gap=$gap: ") {
        fixture.rescuer.rescue(record) shouldBe None
        new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minAnchor=19)
          .rescue(record).isDefined shouldBe true
        record.asSam.getSAMString shouldBe before
      }
    }
  }

  it should "not count ambiguous bases as matches or anchor evidence" in {
    val fixture = new Fixture()
    val record = fixture.read(leftLength=25, rightLength=125, sourceLeft=false)
    record.bases = ("N" * 25) + fixture.right(125)
    fixture.rescuer.rescue(record) shouldBe None
    // N/N must not acquire a positive alignment score from equality in the scorer.
    fixture.refs(0) = fixture.refs(0).patch(575, "N" * 25, 25)
    fixture.rescuer.rescue(record) shouldBe None
    val unknown = fixture.read()
    unknown.bases = "N" * 150
    fixture.rescuer.rescue(unknown) shouldBe None
  }

  it should "exclude an internal ambiguous base from an otherwise sufficient junction anchor" in {
    val fixture = new Fixture()
    for (sourceLeft <- IndexedSeq(true, false)) {
      val record = fixture.sourceRead(sourceLeft, clipLength=20)
      val offset = if (sourceLeft) 135 else 10
      record.bases = record.basesString.updated(offset, 'N')
      val permissive = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minAnchor=19)
      val evidence = permissive.rescue(record).value
      evidence.score shouldBe 285
      evidence.referenceScore should be <= 275
      evidence.identity shouldBe (144.0 / 145)
      evidence.queryCoverage shouldBe 1.0
      val before = record.asSam.getSAMString
      fixture.rescuer.rescue(record) shouldBe None
      evidence.leftAnchor shouldBe (if (sourceLeft) 125 else 19)
      evidence.rightAnchor shouldBe (if (sourceLeft) 19 else 125)
      record.asSam.getSAMString shouldBe before
    }
  }

  it should "enforce identity and whole-query coverage independently" in {
    val fixture = new Fixture()
    val record = fixture.read(leftLength=75, rightLength=75)
    record.bases = record.basesString.updated(35, 'N')
    fixture.rescuer.rescue(record).isDefined shouldBe true
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minIdentity=1.0)
      .rescue(record) shouldBe None
    val clipped = fixture.read(leftLength=75, rightLength=75)
    clipped.bases = ("N" * 16) + clipped.basesString.drop(16)
    fixture.rescuer.rescue(clipped) shouldBe None
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, minQueryCoverage=0.8)
      .rescue(clipped).value.queryCoverage should be < 0.9
  }

  it should "reject equal best query orientations" in {
    val fixture = new Fixture()
    val half = randomBases(75, 91)
    val palindrome = half + reverseComplement(half)
    fixture.refs(0) = fixture.refs(0).patch(525, half, 75)
    fixture.refs(1) = fixture.refs(1).patch(1399, reverseComplement(half), 75)
    val record = fixture.read(leftLength=75, rightLength=75)
    record.basesString shouldBe palindrome
    val aligner = fixtureAligner
    val junction = (fixture.left(200) + fixture.right(200)).getBytes("US-ASCII")
    for (query <- IndexedSeq(palindrome, reverseComplement(palindrome))) {
      val bases = query.getBytes("US-ASCII")
      aligner.align(bases, junction).score shouldBe 300
      aligner.align(bases, junction, 300).size shouldBe 1
      for ((index, cut) <- IndexedSeq((0, 600), (1, 1400))) {
        aligner.align(bases, fixture.reference(index, cut - 205, cut + 205)).score should be <= 290
      }
    }
    fixture.rescuer.rescue(record) shouldBe None
  }

  it should "reject equally optimal junction-spanning endpoint placements" in {
    val fixture = new Fixture()
    val repeat = "ACGTC" * 6
    fixture.refs(0) = fixture.refs(0).patch(570, repeat, repeat.length)
    fixture.refs(1) = fixture.refs(1).patch(1399, repeat, repeat.length)
    val record = fixture.read(leftLength=25, rightLength=25)
    record.basesString shouldBe "ACGTC" * 10
    val aligner = fixtureAligner
    val query = record.bases
    val junction = (fixture.left(200) + fixture.right(200)).getBytes("US-ASCII")
    aligner.align(query, junction).score shouldBe 100
    aligner.align(reverseComplement(record.basesString).getBytes("US-ASCII"), junction).score should be < 100
    val optima = aligner.align(query, junction, 100)
    optima.map(a => (a.queryStart, a.queryEnd, a.targetStart, a.targetEnd)).toSet.size shouldBe 3
    optima.foreach { alignment =>
      alignment.queryStart shouldBe 1
      alignment.queryEnd shouldBe 50
      alignment.cigar.toString shouldBe "50="
      val leftAnchor = 201 - alignment.targetStart
      val rightAnchor = alignment.targetEnd - 200
      leftAnchor should be >= 20
      rightAnchor should be >= 20
    }
    for (bases <- IndexedSeq(query, reverseComplement(record.basesString).getBytes("US-ASCII"));
         (index, cut) <- IndexedSeq((0, 600), (1, 1400))) {
      aligner.align(bases, fixture.reference(index, cut - 205, cut + 205)).score should be <= 62
    }
    fixture.rescuer.rescue(record) shouldBe None
  }

  it should "clamp contig-edge flanks and use the actual incoming length as the boundary" in {
    val fixture = new Fixture(Breakpoint(0, 25, true, 1, 2476, true))
    val rescuer = fixture.rescuer
    rescuer.junctionBoundary shouldBe 25
    val result = rescuer.rescue(fixture.read(leftLength=25, rightLength=25)).value
    result.leftAnchor shouldBe 25
    result.rightAnchor shouldBe 25
    result.targetStart shouldBe 1
    result.targetEnd shouldBe 50
    for (bp <- IndexedSeq(Breakpoint(0, 2476, false, 1, 25, false),
      Breakpoint(0, 25, true, 1, 25, false), Breakpoint(0, 2476, false, 1, 2476, true))) {
      val edge = new Fixture(bp)
      edge.rescuer.junctionBoundary shouldBe 25
      edge.rescuer.rescue(edge.read(leftLength=25, rightLength=25)).value.score shouldBe 100
    }
    val tooShort = new Fixture(Breakpoint(0, 19, true, 1, 1400, true))
    tooShort.rescuer.rescue(tooShort.read(leftLength=19, rightLength=125, sourceLeft=false)) shouldBe None
  }

  it should "cache targets and request complete slop-expanded physical reference neighborhoods" in {
    val fixture = new Fixture()
    val intervals = scala.collection.mutable.ArrayBuffer.empty[(Int, Int, Int)]
    val callback: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
      intervals += ((index, start, end))
      fixture.reference(index, start, end)
    }
    val rescuer = new JunctionRescuer(fixture.bp, callback, fixture.lengths, flank=200)
    rescuer.junctionBoundary shouldBe 200
    intervals.toIndexedSeq shouldBe IndexedSeq((0, 401, 600), (1, 1400, 1599), (0, 395, 805), (1, 1195, 1605))
    rescuer.rescue(fixture.read()).isDefined shouldBe true
    rescuer.rescue(fixture.read()).isDefined shouldBe true
    intervals.size shouldBe 4
    val targets = JunctionRescuer.buildTargets(fixture.bp, fixture.reference, fixture.lengths, flank=200)
    targets.boundary shouldBe 200
    new String(targets.junction, "US-ASCII") shouldBe fixture.left(200) + fixture.right(200)
  }

  it should "uppercase copied query and reference bases without modifying inputs" in {
    for (leftPositive <- IndexedSeq(true, false); rightPositive <- IndexedSeq(true, false)) {
      val fixture = new Fixture(Breakpoint(0, 600, leftPositive, 1, 1400, rightPositive))
      val record = fixture.read()
      val expected = fixture.rescuer.rescue(record).value
      record.bases = record.basesString.toLowerCase(java.util.Locale.ROOT)
      val before = record.asSam.getSAMString
      val fetched = scala.collection.mutable.ArrayBuffer.empty[(Array[Byte], Array[Byte])]
      val reference: (Int, Int, Int) => Array[Byte] = (index, start, end) => {
        val bases = fixture.reference(index, start, end).map(b => (b + 32).toByte)
        fetched += ((bases, bases.clone()))
        bases
      }
      new JunctionRescuer(fixture.bp, reference, fixture.lengths, flank=200, minAnchor=20)
        .rescue(record).value shouldBe expected
      record.asSam.getSAMString shouldBe before
      fetched.foreach { case (bases, original) => bases.toIndexedSeq shouldBe original.toIndexedSeq }
    }
  }

  it should "replace every noncanonical reverse-complemented target base with N" in {
    val bases = "acgtnryswkmbdhv=".getBytes("US-ASCII")
    val before = bases.clone()
    val bp = Breakpoint(0, 1, false, 1, bases.length, false)
    val targets = JunctionRescuer.buildTargets(bp, (_, _, _) => bases, IndexedSeq.fill(2)(bases.length),
      flank=bases.length, slop=0)
    new String(targets.junction, "US-ASCII") shouldBe ("N" * 12 + "ACGT") * 2
    targets.references.foreach { reference =>
      new String(reference, "US-ASCII") shouldBe "ACGTNRYSWKMBDHV="
    }
    bases.toIndexedSeq shouldBe before.toIndexedSeq
  }

  it should "require flank to cover eligible read length plus slop without rejecting irrelevant reads" in {
    val fixture = new Fixture()
    val record = fixture.read()
    val requiredFlank = record.bases.length + 5
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=requiredFlank)
      .rescue(record).value.score shouldBe 300
    val undersized = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=requiredFlank - 1)
    val error = intercept[IllegalArgumentException](undersized.rescue(record))
    error.getMessage should include ("increase --flank")
    error.getMessage should include (s"at least $requiredFlank")
    new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=150, slop=0)
      .rescue(record).value.score shouldBe 300

    val changes: IndexedSeq[SamRecord => Unit] = IndexedSeq(
      r => r.duplicate = true,
      r => r.secondary = true,
      r => r.supplementary = true,
      r => r.pf = false,
      r => r.unmapped = true,
      r => r.mapq = 0,
      r => r.cigar = Cigar("150M"),
      r => r.start += 6)
    changes.foreach { change =>
      val irrelevant = fixture.read()
      change(irrelevant)
      undersized.rescue(irrelevant) shouldBe None
    }

    val hugeSlop = new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, flank=200, slop=Int.MaxValue)
    val overflow = intercept[IllegalArgumentException](hugeSlop.rescue(record))
    overflow.getMessage should include (s"at least ${record.bases.length.toLong + Int.MaxValue.toLong}")
  }

  it should "validate coordinates, configuration, reference lengths and read-length bounds" in {
    val fixture = new Fixture()
    val invalid: IndexedSeq[() => JunctionRescuer] = IndexedSeq(
      () => new JunctionRescuer(fixture.bp.copy(leftRefIndex= -1), fixture.reference, fixture.lengths, 200),
      () => new JunctionRescuer(fixture.bp.copy(rightRefIndex=2), fixture.reference, fixture.lengths, 200),
      () => new JunctionRescuer(fixture.bp.copy(leftPos=0), fixture.reference, fixture.lengths, 200),
      () => new JunctionRescuer(fixture.bp.copy(rightPos=2501), fixture.reference, fixture.lengths, 200),
      () => new JunctionRescuer(fixture.bp, fixture.reference, IndexedSeq(0, 2500), 200),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 0),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 9),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minAnchor=0),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, slop= -1),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minMapq= -1),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minMapq=256),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minScoreDelta= -1),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minIdentity=Double.NaN),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minIdentity=0),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minIdentity=1.1),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minQueryCoverage=Double.NaN),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minQueryCoverage=0),
      () => new JunctionRescuer(fixture.bp, fixture.reference, fixture.lengths, 200, minQueryCoverage=1.1),
      () => new JunctionRescuer(fixture.bp, (_, _, _) => Array.emptyByteArray, fixture.lengths, 200))
    invalid.foreach(create => intercept[IllegalArgumentException](create()))
    val longRead = fixture.read(leftLength=125, rightLength=100)
    intercept[IllegalArgumentException](fixture.rescuer.rescue(longRead))
  }
}
