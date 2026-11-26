"""Polishing's homopolymer stage: every run of the consensus at the length a homopolymer model finds most likely given
the reads' lengths there, each read on its own strand (`_runs`, shasta2's msa1 repair); without a model, ties
between run lengths taken longer (`_lengthen_ties`)."""

from karan import Aligner

from .bands import PolishScratch, TIE_RUN, _Paths, _Realigned, _code
from .homopolymer import FLANK_UNKNOWN, STRAND_UNKNOWN, HomopolymerModel, reads_say_one


def _run_lengths(
    paths: _Paths,
    sequences: List[List[UInt8]],
    members: List[Int],
    strands: List[Int],
    s: Int,
    e: Int,
    base: Int,
    mut lengths: List[Int],
    mut on: List[Int],
    mut which: List[Int],
):
    """The length each member reports for the consensus's run `s ..< e` of `base` (code), its strand, and which
    member it is (`which`, its index in `members`)."""
    var length = paths.stride - 1
    var need_start = s - 1 if s > 0 else 0
    var need_end = e + 1 if e < length else length
    lengths.clear()
    on.clear()
    which.clear()
    for k in range(len(members)):
        if paths.span_start[k] > need_start or paths.span_end[k] < need_end:
            continue
        ref seq = sequences[members[k]]
        var kb = k * paths.stride
        var a = Int(paths.enter[kb + s])
        var b = Int(paths.leave[kb + e])
        while a < b and _code(seq[a]) != base:
            a += 1
        while b > a and _code(seq[b - 1]) != base:
            b -= 1
        var pure = True
        for j in range(a, b):
            pure = pure and _code(seq[j]) == base
        if not pure:
            continue
        if s > 0:
            var at = Int(paths.leave[kb + s - 1])
            if at >= Int(paths.enter[kb + s]) or _code(seq[at]) == base:
                continue
        if e < length:
            var at = Int(paths.leave[kb + e])
            if at >= Int(paths.enter[kb + e + 1]) or _code(seq[at]) == base:
                continue
        lengths.append(b - a)
        which.append(k)
        var r = members[k]
        on.append(strands[r] if r < len(strands) else STRAND_UNKNOWN)


def _lengthen_ties(
    mut aligner: Aligner,
    var current: List[UInt8],
    sequences: List[List[UInt8]],
    unique: List[Int],
    anchors: List[Int],
    weight: List[Int],
    mut realigned: _Realigned,
    mut scratch: PolishScratch,
    threads: Int = 1,
) raises -> List[UInt8]:
    """Every run of one letter at least `TIE_RUN` long one letter longer where that leaves the sequences' summed
    edit distance as it is."""
    ref paths = scratch.paths
    paths.fill(aligner, realigned, current, sequences, unique, anchors, weight, threads=threads)
    var base = paths.total
    ref lengths = scratch.lengths
    ref on = scratch.on
    ref which = scratch.which
    var ties = List[Tuple[Int, Int]]()
    var i = 0
    while i < len(current):
        var e = i + 1
        while e < len(current) and current[e] == current[i]:
            e += 1
        var code = _code(current[i])
        if e - i >= TIE_RUN and code < 4:
            _run_lengths(paths, sequences, unique, List[Int](), i, e, code, lengths, on, which)
            var balance = 0
            for k in range(len(lengths)):
                var w = weight[which[k]] if len(weight) > 0 else 1
                balance += w if lengths[k] > e - i else -w
            if balance == 0 and len(lengths) > 0:
                ties.append((i, e))
        i = e
    var shift = 0
    for t in ties:
        var start = t[0] + shift
        var end = t[1] + shift
        var longer = List[UInt8](capacity=len(current) + 1)
        for k in range(end):
            longer.append(current[k])
        longer.append(current[start])
        for k in range(end, len(current)):
            longer.append(current[k])
        paths.fill(aligner, realigned, longer, sequences, unique, anchors, weight, threads=threads)
        if paths.total == base:
            current = longer^
            shift += 1
    return current^


def _runs(
    mut aligner: Aligner,
    var consensus: List[UInt8],
    sequences: List[List[UInt8]],
    members: List[Int],
    anchors: List[Int],
    model: HomopolymerModel,
    strands: List[Int],
    mut realigned: _Realigned,
    mut scratch: PolishScratch,
    threads: Int = 1,
) raises -> List[UInt8]:
    """Every run of two letters or more in `consensus` at its most likely length under `model`."""
    var length = len(consensus)
    ref paths = scratch.paths
    paths.fill(aligner, realigned, consensus, sequences, members, anchors, threads=threads)
    var out = List[UInt8](capacity=length + length // 8)
    ref lengths = scratch.lengths
    ref on = scratch.on
    ref which = scratch.which
    var s = 0
    while s < length:
        var e = s + 1
        while e < length and consensus[e] == consensus[s]:
            e += 1
        var base = _code(consensus[s])
        if e - s < 2 or base > 3:
            for i in range(s, e):
                out.append(consensus[i])
            s = e
            continue
        var left = _code(consensus[s - 1]) if s > 0 else FLANK_UNKNOWN
        var right = _code(consensus[e]) if e < length else FLANK_UNKNOWN
        _run_lengths(paths, sequences, members, strands, s, e, base, lengths, on, which)
        var n = e - s
        if reads_say_one(lengths):
            n = 1
        elif len(lengths) > 0 and model.covers(base, _median(lengths)):
            var best = model.most_likely(base, left, right, lengths, on)
            if best >= 2:
                n = best
        for _ in range(n):
            out.append(consensus[s])
        s = e
    return out^


def _median(values: List[Int]) -> Int:
    """The lower median of `values` (not empty)."""
    var sorted_values = values.copy()
    sort(sorted_values)
    return sorted_values[(len(sorted_values) - 1) // 2]


def _groups(
    sequences: List[List[UInt8]],
    members: List[Int],
    anchors: List[Int],
    mut unique: List[Int],
    mut weight: List[Int],
    mut same_as: List[Int],
) raises:
    """`members` grouped by letters and anchor: the first of each group in `unique` (in the order they come), how
    many it stands for in `weight`, and per sequence the first of its group in `same_as` (itself if not a member)."""
    var at = Dict[String, Int]()
    same_as = List[Int](capacity=len(sequences))
    for r in range(len(sequences)):
        same_as.append(r)
    for r in members:
        var anchor = anchors[r] if len(anchors) > 0 else 0
        var key = String(anchor) + ":" + String(unsafe_from_utf8=sequences[r])
        var k = at.get(key, -1)
        if k < 0:
            at[key] = len(unique)
            unique.append(r)
            weight.append(1)
        else:
            weight[k] += 1
            same_as[r] = unique[k]
