"""The backbone of a pairwise alignment: chained anchors, refined inside the gaps they leave.

A level finds anchors in a stretch of both sequences, chains them between the stretch's fixed ends and cuts the
overlaps away. Every gap between them still longer than the budget is searched again, its anchors found and
counted inside the gap alone, where a string repeated across the whole sequence may be unique (LAGAN's
recursion): at the next level, or with `stay` at the same level until a gap yields nothing, and only then at the
next. What is left between the anchors are the windows an exact aligner fills.
"""

from std.collections import Dict

from .alphabet import encode_base
from .anchors import Anchor, find_anchors, mums, mums_or_smems, rare, rare_or_maximal, smems
from .chain import ChainCosts, chain, without_overlaps


comptime MUMS = 0
"""Anchors whose string occurs once in each sequence."""
comptime REFERENCE_UNIQUE = 1
"""Anchors whose string occurs once in the partner, any number of times in the query (MUMmer's -mumreference)."""
comptime RARE = 2
"""Anchors whose string occurs at most `Level.cap` times in each."""
comptime SMEMS = 3
"""Anchors no other anchor's query interval strictly contains."""
comptime ALL = 4
"""Every MEM."""
comptime RARE_OR_MAXIMAL = 5
"""Anchors that are rare (`Level.cap` copies at most, MUMs among them) or SMEMs."""
comptime MUMS_OR_SMEMS = 6
"""MUMs and SMEMs, no other MEM."""


struct Level(ImplicitlyCopyable):
    """One round of anchoring: which anchors, of how many letters at least, counted to which cap."""

    var kind: Int
    var min_len: Int
    var cap: Int

    def __init__(out self, kind: Int, min_len: Int, cap: Int = 1):
        self.kind = kind
        self.min_len = min_len
        self.cap = cap

    def select(self, var anchors: List[Anchor]) -> List[Anchor]:
        if self.kind == MUMS:
            return mums(anchors)
        if self.kind == REFERENCE_UNIQUE:
            var out = List[Anchor]()
            for a in anchors:
                if a.ref_copies == 1:
                    out.append(a)
            return out^
        if self.kind == RARE:
            return rare(anchors, self.cap)
        if self.kind == SMEMS:
            return smems(anchors)
        if self.kind == RARE_OR_MAXIMAL:
            return rare_or_maximal(anchors, self.cap)
        if self.kind == MUMS_OR_SMEMS:
            return mums_or_smems(anchors)
        return anchors^

    def counted_to(self) -> Int:
        """The cap its anchors are counted to: past `cap`, so copies above it show (and weigh)."""
        return max(self.cap, 1) + 1 if self.kind == RARE or self.kind == RARE_OR_MAXIMAL else 2


struct Window(ImplicitlyCopyable):
    """The letters between two backbone anchors: `ref_start ..< ref_end` and `query_start ..< query_end`."""

    var ref_start: Int
    var ref_end: Int
    var query_start: Int
    var query_end: Int

    def __init__(out self, ref_start: Int, ref_end: Int, query_start: Int, query_end: Int):
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.query_start = query_start
        self.query_end = query_end

    def size(self) -> Int:
        return max(self.ref_end - self.ref_start, self.query_end - self.query_start)


struct RepeatMask(Movable):
    """The stretches of a sequence inside tandem repeats: a k-mer start is masked when the same k-mer starts
    again within `distance` letters. An anchor any of whose k-mers is masked sits in a repeat, where an exact
    match, unique or not, may pair different copies."""

    var k: Int
    var masked_before: List[Int32]
    """How many masked k-mer starts precede each position."""

    def __init__(out self, sequence: ImmSpan[UInt8, _], k: Int, distance: Int):
        self.k = k
        var n = len(sequence)
        self.masked_before = List[Int32](length=n + 1, fill=0)
        if k <= 0 or n < k:
            return
        var marked = List[Bool](length=n, fill=False)
        var last = Dict[UInt64, Int]()
        var mask = (UInt64(1) << UInt64(2 * k)) - 1
        var key = UInt64(0)
        var valid = 0
        for i in range(n):
            var code = encode_base(sequence[i])
            valid = valid + 1 if code < 4 else 0
            key = ((key << 2) | UInt64(code & 3)) & mask
            if valid < k:
                continue
            var start = i - k + 1
            var seen = last.get(key)
            if seen and start - seen.value() <= distance:
                marked[start] = True
                marked[seen.value()] = True
            last[key] = start
        for i in range(n):
            self.masked_before[i + 1] = self.masked_before[i] + Int32(1 if marked[i] else 0)

    def covers(self, start: Int, end: Int) -> Bool:
        """Whether any k-mer of `start ..< end` is masked."""
        if self.k <= 0 or end - start < self.k:
            return False
        return self.masked_before[end - self.k + 1] > self.masked_before[start]


def backbone(
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    levels: List[Level],
    costs: ChainCosts,
    budget: Int,
    stay: Bool = False,
    max_anchors: Int = 50_000,
    lookback: Int = 0,
    repeat_k: Int = 0,
    repeat_distance: Int = 0,
    mask_depth: Int = Int.MAX,
) raises -> List[Anchor]:
    """Non-overlapping anchors, in order on both sequences: level 0 over the whole pair, then inside every gap
    still longer than `budget` letters the next level, or with `stay` the same one until a gap yields nothing.
    A gap with more than `max_anchors` anchors is left to the next level; `lookback` caps how far back each
    anchor's chain search goes (0: no cap). With `repeat_k`, the anchors inside tandem repeats of the partner
    (`RepeatMask` of `repeat_k`-mers within `repeat_distance`) are dropped before chaining, at the first
    `mask_depth` depths of the recursion (every depth by default): with 1, the whole pair's anchors keep out of
    the repeats and so pin the flanks around them, and inside a gap between anchors, a repeat array, a match is
    unique within the gap's two pieces, as a copy that carries something of its own is."""
    var out = List[Anchor]()
    var repeats = RepeatMask(partner, repeat_k, repeat_distance)
    _refine(partner, query, 0, 0, levels, 0, costs, budget, stay, max_anchors, lookback, repeats, mask_depth, 0, out)
    return out^


def _refine(
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    ref_offset: Int,
    query_offset: Int,
    levels: List[Level],
    level: Int,
    costs: ChainCosts,
    budget: Int,
    stay: Bool,
    max_anchors: Int,
    lookback: Int,
    repeats: RepeatMask,
    mask_depth: Int,
    depth: Int,
    mut out: List[Anchor],
) raises:
    if level >= len(levels) or max(len(partner), len(query)) <= budget:
        return
    var current = levels[level]
    var found = current.select(find_anchors(partner, query, current.min_len, current.counted_to()))
    if repeats.k > 0 and depth < mask_depth:
        var outside = List[Anchor]()
        for a in found:
            if not repeats.covers(ref_offset + a.ref_pos, ref_offset + a.ref_end()):
                outside.append(a)
        found = outside^
    var picked = List[Anchor]()
    if len(found) <= max_anchors:
        picked = without_overlaps(chain(found, len(partner), len(query), costs, lookback).anchors)
    var next = level if stay else level + 1
    if len(picked) == 0:
        _refine(
            partner,
            query,
            ref_offset,
            query_offset,
            levels,
            level + 1,
            costs,
            budget,
            stay,
            max_anchors,
            lookback,
            repeats,
            mask_depth,
            depth,
            out,
        )
        return
    var ref_at = 0
    var query_at = 0
    for a in picked:
        _refine(
            partner[ref_at : a.ref_pos],
            query[query_at : a.query_pos],
            ref_offset + ref_at,
            query_offset + query_at,
            levels,
            next,
            costs,
            budget,
            stay,
            max_anchors,
            lookback,
            repeats,
            mask_depth,
            depth + 1,
            out,
        )
        out.append(a.shifted(ref_offset, query_offset))
        ref_at = a.ref_end()
        query_at = a.query_end()
    _refine(
        partner[ref_at:],
        query[query_at:],
        ref_offset + ref_at,
        query_offset + query_at,
        levels,
        next,
        costs,
        budget,
        stay,
        max_anchors,
        lookback,
        repeats,
        mask_depth,
        depth + 1,
        out,
    )


def windows(anchors: List[Anchor], ref_len: Int, query_len: Int) -> List[Window]:
    """The windows between non-overlapping `anchors` in order, from the start of both sequences to their ends,
    empty ones included."""
    var out = List[Window]()
    var ref_at = 0
    var query_at = 0
    for a in anchors:
        out.append(Window(ref_at, a.ref_pos, query_at, a.query_pos))
        ref_at = a.ref_end()
        query_at = a.query_end()
    out.append(Window(ref_at, ref_len, query_at, query_len))
    return out^
