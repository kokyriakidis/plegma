"""CIGAR parsing and alignment coordinate walkers.

Direct byte-level CIGAR scanning avoids allocating intermediate `(op, count)` lists in hot loops while keeping the
public `Runs` / `cigar_runs` representation for callers that need it.
"""

from karan import Alignment

comptime MATCH = UInt8(ord("="))
comptime MISMATCH = UInt8(ord("X"))
comptime INSERTION = UInt8(ord("I"))
comptime DELETION = UInt8(ord("D"))
comptime GAP = UInt8(ord("-"))
comptime ALIGN_MATCH = UInt8(ord("M"))


struct Runs(Movable):
    """A CIGAR's runs: their letters and lengths."""

    var ops: List[UInt8]
    var counts: List[Int]

    def __init__(out self):
        self.ops = List[UInt8]()
        self.counts = List[Int]()


def cigar_runs(cigar: StringSlice) raises -> Runs:
    """`cigar`'s runs."""
    var runs = Runs()
    var length = 0
    for c in cigar.as_bytes():
        if c >= 48 and c <= 57:
            length = length * 10 + Int(c - 48)
            continue
        if length == 0:
            raise Error(t"run without a length in {cigar}")
        runs.ops.append(c)
        runs.counts.append(length)
        length = 0
    return runs^


def pairing_from_cigar(cigar: StringSlice, q_length: Int, ref_offset: Int = 0) -> List[Int]:
    """For each query letter of a global alignment's CIGAR, the reference letter it pairs with (-1: inserted)."""
    var out = List[Int](capacity=q_length)
    append_pairing_from_cigar(cigar, ref_offset, out)
    return out^


def append_pairing_from_cigar(cigar: StringSlice, ref_offset: Int, mut out: List[Int]):
    """Appends query-to-reference pairings from `cigar` (starting at `ref_offset`) into `out`."""
    var i = ref_offset
    var run_n = 0
    for ch in cigar.as_bytes():
        if ch >= 48 and ch <= 57:
            run_n = run_n * 10 + Int(ch - 48)
            continue
        var n = run_n if run_n > 0 else 1
        run_n = 0
        if ch == DELETION:
            i += n
        elif ch == INSERTION:
            for _ in range(n):
                out.append(-1)
        else:
            for _ in range(n):
                out.append(i)
                i += 1


def query_at(alignment: Alignment) -> List[Int]:
    """For each reference letter `reference_start ..= reference_end` of `alignment`, the query letter it stands
    before (a gap's letters count before the reference letter after them): the query's letters over reference
    `s ..< e` are `query_at[s - reference_start] ..< query_at[e - reference_start]`."""
    var out = List[Int](capacity=alignment.reference_end - alignment.reference_start + 1)
    var q = alignment.query_start
    var run_n = 0
    for ch in alignment.cigar.as_bytes():
        if ch >= 48 and ch <= 57:
            run_n = run_n * 10 + Int(ch - 48)
            continue
        var n = run_n if run_n > 0 else 1
        run_n = 0
        if ch == INSERTION:
            q += n
            continue
        if ch == DELETION:
            for _ in range(n):
                out.append(q)
        else:
            for _ in range(n):
                out.append(q)
                q += 1
    out.append(q)
    return out^


def haplotype_at(alignment: Alignment) -> List[Int]:
    """`query_at` the other way round: for a haplotype placed inside a read (`alignment`'s query the haplotype, its
    reference the read; `Mode.INFIX`, so every haplotype letter is placed), for each haplotype letter `0 ..= len` the
    read letter it stands before, a gap's letters counting before the haplotype letter after them."""
    var out = List[Int](capacity=alignment.query_end - alignment.query_start + 1)
    var r = alignment.reference_start
    var run_n = 0
    for ch in alignment.cigar.as_bytes():
        if ch >= 48 and ch <= 57:
            run_n = run_n * 10 + Int(ch - 48)
            continue
        var n = run_n if run_n > 0 else 1
        run_n = 0
        if ch == DELETION:
            r += n
            continue
        if ch == INSERTION:
            for _ in range(n):
                out.append(r)
        else:
            for _ in range(n):
                out.append(r)
                r += 1
    out.append(r)
    return out^
