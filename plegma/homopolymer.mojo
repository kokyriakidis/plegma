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

    def _done(mut self, bytes: Span[UInt8, _]) raises:
        """Past the field's end: a comma, or the line's end on the last field."""
        if self.at < self.end:
            if bytes[self.at] != UInt8(ord(",")):
                raise Error("homopolymer model: a malformed field")
            self.at += 1

    def flank(mut self, bytes: Span[UInt8, _]) raises -> Int:
        if self.at >= self.end:
            raise Error("homopolymer model: a line short of fields")
        var got = _flank(bytes[self.at])
        self.at += 1
        self._done(bytes)
        return got

    def integer(mut self, bytes: Span[UInt8, _]) raises -> Int:
        var got = 0
        var digits = 0
        while self.at < self.end and bytes[self.at] >= UInt8(ord("0")) and bytes[self.at] <= UInt8(ord("9")):
            got = got * 10 + Int(bytes[self.at] - UInt8(ord("0")))
            self.at += 1
            digits += 1
        if digits == 0:
            raise Error("homopolymer model: a number expected")
        self._done(bytes)
        return got

    def skip(mut self, bytes: Span[UInt8, _]) raises:
        while self.at < self.end and bytes[self.at] != UInt8(ord(",")):
            self.at += 1
        self._done(bytes)

    def decimal(mut self, bytes: Span[UInt8, _]) raises -> Float64:
        """A decimal number: a sign, digits, a fraction and an exponent, each optional but the digits."""
        var negative = False
        if self.at < self.end and (bytes[self.at] == UInt8(ord("-")) or bytes[self.at] == UInt8(ord("+"))):
            negative = bytes[self.at] == UInt8(ord("-"))
            self.at += 1
        var value = 0.0
        var digits = 0
        while self.at < self.end and bytes[self.at] >= UInt8(ord("0")) and bytes[self.at] <= UInt8(ord("9")):
            value = value * 10.0 + Float64(Int(bytes[self.at] - UInt8(ord("0"))))
            self.at += 1
            digits += 1
        if self.at < self.end and bytes[self.at] == UInt8(ord(".")):
            self.at += 1
            var scale = 0.1
            while self.at < self.end and bytes[self.at] >= UInt8(ord("0")) and bytes[self.at] <= UInt8(ord("9")):
                value += scale * Float64(Int(bytes[self.at] - UInt8(ord("0"))))
                scale *= 0.1
                self.at += 1
                digits += 1
        if digits == 0:
            raise Error("homopolymer model: a number expected")
        if self.at < self.end and (bytes[self.at] == UInt8(ord("e")) or bytes[self.at] == UInt8(ord("E"))):
            self.at += 1
            var minus = self.at < self.end and bytes[self.at] == UInt8(ord("-"))
            if self.at < self.end and (bytes[self.at] == UInt8(ord("-")) or bytes[self.at] == UInt8(ord("+"))):
                self.at += 1
            var power = 0
            while self.at < self.end and bytes[self.at] >= UInt8(ord("0")) and bytes[self.at] <= UInt8(ord("9")):
                power = power * 10 + Int(bytes[self.at] - UInt8(ord("0")))
                self.at += 1
            value *= 10.0 ** Float64(-power if minus else power)
        self._done(bytes)
        return -value if negative else value


struct HomopolymerModel(Copyable, Movable):
    """A homopolymer model read from shasta2's csv (`load`); empty (`loaded` false) by default. Log probabilities in
    dB, as the file writes them."""

    var n_count: Int
    var m_count: Int
    var table: List[Float64]
    """log P(m | n, ...) at `_at(base, strand, left, right, n) + m`."""
    var present: List[Bool]
    """Whether the row of (base, strand, left, right, n) is in the model, at `_row(...)`."""
    var log_prior: List[Float64]

    def __init__(out self):
        self.n_count = 0
        self.m_count = 0
        self.table = List[Float64]()
        self.present = List[Bool]()
        self.log_prior = List[Float64]()

    @staticmethod
    def load(path: String) raises -> HomopolymerModel:
        """The model in shasta2's csv at `path`; every row checked whole, as shasta2 checks it."""
        var text: String
        with open(path, "r") as f:
            text = f.read()
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
            var base = fields.flank(bytes)
            var strand = fields.integer(bytes)
            var left = fields.flank(bytes)
            var right = fields.flank(bytes)
            var n = fields.integer(bytes)
            var m = fields.integer(bytes)
            var count = fields.integer(bytes)
            fields.skip(bytes)
            var value = fields.decimal(bytes)
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
        return model^

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

    def most_likely(self, base: Int, left: Int, right: Int, lengths: List[Int], strands: List[Int]) -> Int:
        """The run length n most likely given the reads' lengths `lengths` (each on `strands`' strand: 0, 1 or
        `STRAND_UNKNOWN`), `base` between `left` and `right` (0..3, or `FLANK_UNKNOWN`): the n that maximizes
        log P(n) plus the sum of log P(m | n, ...), the shorter on ties; -1 when no n is in the model."""
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
