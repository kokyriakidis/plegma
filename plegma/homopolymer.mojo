"""An empirical model of homopolymer length errors, as shasta2 trains one: the probability P(m | n, base, strand,
left, right) that a read reports a run of `base` m long where it is n long, between the letters `left` and `right`,
the read on `strand`. With it a run's most likely length given the lengths the reads report replaces their vote,
which nanopore reads bias short (shasta2's msa1 repair: on its chr21 assembly two thirds of the errors the plain
consensus makes).

The file is shasta2's csv (conf/HomopolymerModel-*.csv in github.com/paoloshasta/shasta2):

    base,strand,left,right,n,m,count,probability,logProbabilityDb

one line per (base, strand, left, right, n, m); left and right are A, C, G, T or `*` (the row for any flanks, required:
used where the flanks have no row of their own); n from 2; logProbabilityDb is 10 log10 of the probability; the
largest n and m each stand for that length or more. The prior P(n) is the `*` rows' counts per n, pooled over bases
and strands, one added to each n present. Strand 0 is a read as the consensus runs, 1 one reverse-complemented to it;
a read of unknown strand (`STRAND_UNKNOWN`) weighs the two strands' probabilities equally. As in shasta2, an n is
possible where both strands have its row.
"""

from std.math import log10
from std.sys import simd_width_of

comptime HEADER = "base,strand,left,right,n,m,count,probability,logProbabilityDb"
comptime STRAND_UNKNOWN = 2
"""A read's strand when its input does not say."""
comptime FLANK_UNKNOWN = 4
"""A run's flank at the consensus's end, or not a letter of ACGT."""
comptime _STRANDS = 2
comptime _FLANKS = 5


def _flank(code: UInt8) raises -> Int:
    """A, C, G, T as 0..3, `*` as `FLANK_UNKNOWN`."""
    if code == UInt8(ord("A")):
        return 0
    if code == UInt8(ord("C")):
        return 1
    if code == UInt8(ord("G")):
        return 2
    if code == UInt8(ord("T")):
        return 3
    if code == UInt8(ord("*")):
        return FLANK_UNKNOWN
    raise Error("homopolymer model: a flank or base not A, C, G, T or *")


struct _Fields:
    """The csv's fields read in place, byte by byte (the model is 200,000 lines; a caller runs plegma once per local
    assembly, so loading it must take milliseconds)."""

    var at: Int
    var end: Int

    def __init__(out self, at: Int, end: Int):
        self.at = at
        self.end = end

    def _done(mut self, bytes: Pointer[UInt8, _]) raises:
        """Past the field's end: a comma, or the line's end on the last field."""
        if self.at < self.end:
            if bytes[unsafe_offset=self.at] != UInt8(ord(",")):
                raise Error("homopolymer model: a malformed field")
            self.at += 1

    def flank(mut self, bytes: Pointer[UInt8, _]) raises -> Int:
        if self.at >= self.end:
            raise Error("homopolymer model: a line short of fields")
        var got = _flank(bytes[unsafe_offset=self.at])
        self.at += 1
        self._done(bytes)
        return got

    def integer(mut self, bytes: Pointer[UInt8, _]) raises -> Int:
        var got = 0
        var digits = 0
        while (
            self.at < self.end
            and bytes[unsafe_offset=self.at] >= UInt8(ord("0"))
            and bytes[unsafe_offset=self.at] <= UInt8(ord("9"))
        ):
            got = got * 10 + Int(bytes[unsafe_offset=self.at] - UInt8(ord("0")))
            self.at += 1
            digits += 1
        if digits == 0:
            raise Error("homopolymer model: a number expected")
        self._done(bytes)
        return got

    def skip(mut self, bytes: Pointer[UInt8, _]) raises:
        while self.at < self.end and bytes[unsafe_offset=self.at] != UInt8(ord(",")):
            self.at += 1
        self._done(bytes)

    def decimal(mut self, bytes: Pointer[UInt8, _]) raises -> Float64:
        """A decimal number: a sign, digits, a fraction and an exponent, each optional but the digits."""
        var negative = False
        if self.at < self.end and (
            bytes[unsafe_offset=self.at] == UInt8(ord("-")) or bytes[unsafe_offset=self.at] == UInt8(ord("+"))
        ):
            negative = bytes[unsafe_offset=self.at] == UInt8(ord("-"))
            self.at += 1
        var value = 0.0
        var digits = 0
        while (
            self.at < self.end
            and bytes[unsafe_offset=self.at] >= UInt8(ord("0"))
            and bytes[unsafe_offset=self.at] <= UInt8(ord("9"))
        ):
            value = value * 10.0 + Float64(Int(bytes[unsafe_offset=self.at] - UInt8(ord("0"))))
            self.at += 1
            digits += 1
        if self.at < self.end and bytes[unsafe_offset=self.at] == UInt8(ord(".")):
            self.at += 1
            var scale = 0.1
            while (
                self.at < self.end
                and bytes[unsafe_offset=self.at] >= UInt8(ord("0"))
                and bytes[unsafe_offset=self.at] <= UInt8(ord("9"))
            ):
                value += scale * Float64(Int(bytes[unsafe_offset=self.at] - UInt8(ord("0"))))
                scale *= 0.1
                self.at += 1
                digits += 1
        if digits == 0:
            raise Error("homopolymer model: a number expected")
        if self.at < self.end and (
            bytes[unsafe_offset=self.at] == UInt8(ord("e")) or bytes[unsafe_offset=self.at] == UInt8(ord("E"))
        ):
            self.at += 1
            var minus = self.at < self.end and bytes[unsafe_offset=self.at] == UInt8(ord("-"))
            if self.at < self.end and (
                bytes[unsafe_offset=self.at] == UInt8(ord("-")) or bytes[unsafe_offset=self.at] == UInt8(ord("+"))
            ):
                self.at += 1
            var power = 0
            while (
                self.at < self.end
                and bytes[unsafe_offset=self.at] >= UInt8(ord("0"))
                and bytes[unsafe_offset=self.at] <= UInt8(ord("9"))
            ):
                power = power * 10 + Int(bytes[unsafe_offset=self.at] - UInt8(ord("0")))
                self.at += 1
            value *= 10.0 ** Float64(-power if minus else power)
        self._done(bytes)
        return -value if negative else value


struct HomopolymerModel(Copyable, Movable):
    """A homopolymer model read from shasta2's csv (`load` a file, `parse` its text); empty (`loaded` false) by default. Log probabilities in
    dB, as the file writes them."""

    var n_count: Int
    var m_count: Int
    var table: List[Float64]
    """log P(m | n, ...) at `_at(base, strand, left, right, n) + m`."""
    var present: List[Bool]
    """Whether the row of (base, strand, left, right, n) is in the model, at `_row(...)`."""
    var log_prior: List[Float64]
    var n_pad: Int
    """`n_count` rounded up to a vector."""
    var by_n: List[Float64]
    """log P(m | n, ...) for every n at once: the row (base, strand, left, right, m) from `_by_n(...)`, n_pad
    long, each n's value as `_offset` resolves it (its flanks' row, else the `*` row; 0 where neither)."""
    var start: List[Float64]
    """Per (base, left, right) and n: log P(n), or -Float64.MAX where n is not possible (not in the model for
    either strand)."""

    def __init__(out self):
        self.n_count = 0
        self.m_count = 0
        self.table = List[Float64]()
        self.present = List[Bool]()
        self.log_prior = List[Float64]()
        self.n_pad = 0
        self.by_n = List[Float64]()
        self.start = List[Float64]()

    @staticmethod
    def load(path: String) raises -> HomopolymerModel:
        """The model in shasta2's csv at `path`; every row checked whole, as shasta2 checks it."""
        var text: String
        with open(path, "r") as f:
            text = f.read()
        return HomopolymerModel.parse(text, path)

    @staticmethod
    def parse(text: String, source: String = "text") raises -> HomopolymerModel:
        """The model in `text`, shasta2's csv (a caller holding it in memory, as shasta2 embeds its models); `source`
        names it in errors."""
        var path = source
        var bytes = text.as_bytes()
        var line_start = 0
        var line_number = 0
        # A line is some 30 bytes.
        var lines = len(bytes) // 24 + 1
        var keys = List[Int](capacity=lines)
        var ms = List[Int](capacity=lines)
        var counts = List[Int](capacity=lines)
        var values = List[Float64](capacity=lines)
        var n_count = 0
        var m_count = 0
        while line_start < len(bytes):
            var line_end = line_start
            while line_end < len(bytes) and bytes[line_end] != UInt8(ord("\n")):
                line_end += 1
            var next = line_end + 1
            if line_end > line_start and bytes[line_end - 1] == UInt8(ord("\r")):
                line_end -= 1
            line_number += 1
            if line_number == 1:
                if String(unsafe_from_utf8=bytes[line_start:line_end]) != HEADER:
                    raise Error(t"homopolymer model {path}: not shasta2's csv header")
                line_start = next
                continue
            if line_end == line_start:
                line_start = next
                continue
            var fields = _Fields(line_start, line_end)
            var at = bytes.unsafe_ptr()
            var base = fields.flank(at)
            var strand = fields.integer(at)
            var left = fields.flank(at)
            var right = fields.flank(at)
            var n = fields.integer(at)
            var m = fields.integer(at)
            var count = fields.integer(at)
            fields.skip(at)
            var value = fields.decimal(at)
            if fields.at != line_end:
                raise Error(t"homopolymer model {path}: line {line_number}: 9 fields expected")
            if base == FLANK_UNKNOWN or strand > 1 or n < 2 or not (value <= 0.0):
                raise Error(t"homopolymer model {path}: line {line_number}: out of range")
            n_count = max(n_count, n + 1)
            m_count = max(m_count, m + 1)
            keys.append((((base * _STRANDS + strand) * _FLANKS + left) * _FLANKS + right) << 16 | n)
            ms.append(m)
            counts.append(count)
            values.append(value)
            line_start = next
        if len(keys) == 0:
            raise Error(t"homopolymer model {path}: empty")
        var model = HomopolymerModel()
        model.n_count = n_count
        model.m_count = m_count
        var rows = 4 * _STRANDS * _FLANKS * _FLANKS * n_count
        model.present = List[Bool](length=rows, fill=False)
        model.table = List[Float64](length=rows * m_count, fill=1.0)
        var prior = List[Float64](length=n_count, fill=0.0)
        for k in range(len(keys)):
            var row = (keys[k] >> 16) * n_count + (keys[k] & 0xFFFF)
            model.present[row] = True
            var at = row * m_count + ms[k]
            if model.table[at] <= 0.0:
                raise Error(t"homopolymer model {path}: a line twice")
            model.table[at] = values[k]
            # Any flanks (`*` on both sides): the prior's counts.
            if (keys[k] >> 16) % (_FLANKS * _FLANKS) == FLANK_UNKNOWN * _FLANKS + FLANK_UNKNOWN:
                prior[keys[k] & 0xFFFF] += Float64(counts[k])
        # Every row present is whole; a flanks' row needs its `*` row; every base and strand has `*` rows.
        for base in range(4):
            for strand in range(2):
                var any_star = False
                for left in range(_FLANKS):
                    for right in range(_FLANKS):
                        for n in range(n_count):
                            if not model.present[model._row(base, strand, left, right, n)]:
                                continue
                            var star = model._row(base, strand, FLANK_UNKNOWN, FLANK_UNKNOWN, n)
                            if left == FLANK_UNKNOWN and right == FLANK_UNKNOWN:
                                any_star = True
                            elif not model.present[star]:
                                raise Error(t"homopolymer model {path}: a flanks' row without its * row")
                            var at = model._at(base, strand, left, right, n)
                            for m in range(m_count):
                                if model.table[at + m] > 0.0:
                                    raise Error(t"homopolymer model {path}: a row missing a line")
                if not any_star:
                    raise Error(t"homopolymer model {path}: a base and strand without * rows")
        var total = 0.0
        for n in range(n_count):
            var any_base = False
            for base in range(4):
                any_base = any_base or model.present[model._row(base, 0, FLANK_UNKNOWN, FLANK_UNKNOWN, n)]
            if any_base:
                prior[n] += 1.0
            total += prior[n]
        model.log_prior = List[Float64](length=n_count, fill=-Float64.MAX)
        for n in range(n_count):
            if prior[n] > 0.0:
                model.log_prior[n] = 10.0 * log10(prior[n] / total)
        model._resolve()
        return model^

    def _resolve(mut self):
        """`by_n` and `start` from `table`, `present` and `log_prior`."""
        comptime W = simd_width_of[DType.float64]()
        self.n_pad = (self.n_count + W - 1) // W * W
        self.by_n = List[Float64](length=4 * 2 * _FLANKS * _FLANKS * self.m_count * self.n_pad, fill=0.0)
        self.start = List[Float64](length=4 * _FLANKS * _FLANKS * self.n_pad, fill=-Float64.MAX)
        for base in range(4):
            for left in range(_FLANKS):
                for right in range(_FLANKS):
                    for n in range(self.n_count):
                        var at0 = self._offset(base, 0, left, right, n)
                        var at1 = self._offset(base, 1, left, right, n)
                        if self.log_prior[n] != -Float64.MAX and at0 >= 0 and at1 >= 0:
                            self.start[((base * _FLANKS + left) * _FLANKS + right) * self.n_pad + n] = self.log_prior[n]
                        for strand in range(2):
                            var at = at0 if strand == 0 else at1
                            if at < 0:
                                continue
                            for m in range(self.m_count):
                                self.by_n[self._by_n(base, strand, left, right, m) + n] = self.table[at + m]

    def _by_n(self, base: Int, strand: Int, left: Int, right: Int, m: Int) -> Int:
        return ((((base * 2 + strand) * _FLANKS + left) * _FLANKS + right) * self.m_count + m) * self.n_pad

    def loaded(self) -> Bool:
        return self.n_count > 0

    def _row(self, base: Int, strand: Int, left: Int, right: Int, n: Int) -> Int:
        return (((base * _STRANDS + strand) * _FLANKS + left) * _FLANKS + right) * self.n_count + n

    def _at(self, base: Int, strand: Int, left: Int, right: Int, n: Int) -> Int:
        return self._row(base, strand, left, right, n) * self.m_count

    def _offset(self, base: Int, strand: Int, left: Int, right: Int, n: Int) -> Int:
        """Where the row for these flanks starts in `table`, else the `*` row's; -1 when n is not in the model."""
        var row = self._row(base, strand, left, right, n)
        if self.present[row]:
            return row * self.m_count
        row = self._row(base, strand, FLANK_UNKNOWN, FLANK_UNKNOWN, n)
        return row * self.m_count if self.present[row] else -1

    def covers(self, base: Int, n: Int) -> Bool:
        """Whether runs of `base` n long are in the model (both strands' `*` rows)."""
        if base < 0 or base > 3 or n < 0 or n >= self.n_count:
            return False
        return (
            self.present[self._row(base, 0, FLANK_UNKNOWN, FLANK_UNKNOWN, n)]
            and self.present[self._row(base, 1, FLANK_UNKNOWN, FLANK_UNKNOWN, n)]
        )

    def log_probability(self, base: Int, left: Int, right: Int, n: Int, m: Int, strand: Int) -> Float64:
        """10 log10 P(m | n): a read on `strand` (0, 1 or `STRAND_UNKNOWN`, the two strands' probabilities
        half each) reporting a run of `base` between `left` and `right` m long where it is n long (the flanks' row,
        else the `*` row); -Float64.MAX when n is not in the model."""
        if n < 0 or n >= self.n_count:
            return -Float64.MAX
        var at0 = self._offset(base, 0, left, right, n)
        var at1 = self._offset(base, 1, left, right, n)
        if at0 < 0 or at1 < 0:
            return -Float64.MAX
        var k = min(max(m, 0), self.m_count - 1)
        if strand == 0:
            return self.table[at0 + k]
        if strand == 1:
            return self.table[at1 + k]
        return 10.0 * log10(0.5 * (10.0 ** (self.table[at0 + k] / 10.0) + 10.0 ** (self.table[at1 + k] / 10.0)))

    def most_likely(self, base: Int, left: Int, right: Int, lengths: List[Int], strands: List[Int]) -> Int:
        """The run length n most likely given the reads' lengths `lengths` (each on `strands`' strand: 0, 1 or
        `STRAND_UNKNOWN`), `base` between `left` and `right` (0..3, or `FLANK_UNKNOWN`): the n that maximizes
        log P(n) plus the sum of log P(m | n, ...), the shorter on ties; -1 when no n is in the model."""
        var known = True
        for strand in strands:
            known = known and (strand == 0 or strand == 1)
        if known:
            # Every n at once, a vector at a time; each n's sum in the reads' order, as below.
            comptime W = simd_width_of[DType.float64]()
            var first = ((base * _FLANKS + left) * _FLANKS + right) * self.n_pad
            var scores = List[Float64](capacity=self.n_pad)
            for n in range(self.n_pad):
                scores.append(self.start[first + n])
            var out = scores.unsafe_ptr()
            var rows = self.by_n.unsafe_ptr()
            for k in range(len(lengths)):
                var row = rows.unsafe_offset(
                    self._by_n(base, strands[k], left, right, min(lengths[k], self.m_count - 1))
                )
                var n = 0
                while n < self.n_pad:
                    out.unsafe_store(n, out.unsafe_load[width=W](n) + row.unsafe_load[width=W](n))
                    n += W
            var best_n = -1
            var best_value = -Float64.MAX
            for n in range(self.n_count):
                if self.start[first + n] == -Float64.MAX:
                    continue
                if scores[n] > best_value:
                    best_n = n
                    best_value = scores[n]
            return best_n
        var best = -1
        var best_score = -Float64.MAX
        for n in range(self.n_count):
            var score = self.log_prior[n]
            if score == -Float64.MAX:
                continue
            var at0 = self._offset(base, 0, left, right, n)
            var at1 = self._offset(base, 1, left, right, n)
            if at0 < 0 or at1 < 0:
                continue
            for k in range(len(lengths)):
                var m = min(lengths[k], self.m_count - 1)
                if strands[k] == 0:
                    score += self.table[at0 + m]
                elif strands[k] == 1:
                    score += self.table[at1 + m]
                else:
                    # Either strand, half each.
                    score += 10.0 * log10(
                        0.5 * (10.0 ** (self.table[at0 + m] / 10.0) + 10.0 ** (self.table[at1 + m] / 10.0))
                    )
            if score > best_score:
                best = n
                best_score = score
        return best
