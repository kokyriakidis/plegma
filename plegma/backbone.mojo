"""The backbone of a pairwise alignment: chained anchors, refined inside the gaps they leave.

A level finds anchors in a stretch of both sequences, chains them between the stretch's fixed ends and cuts the
overlaps away. Every gap between them still longer than the budget is searched again, its anchors found and
counted inside the gap alone, where a string repeated across the whole sequence may be unique (LAGAN's
recursion): at the next level, or with `stay` at the same level until a gap yields nothing, and only then at the
next. What is left between the anchors are the windows an exact aligner fills.
"""

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


def backbone(
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    levels: List[Level],
    costs: ChainCosts,
    budget: Int,
    stay: Bool = False,
    max_anchors: Int = 50_000,
    lookback: Int = 0,
) raises -> List[Anchor]:
    """Non-overlapping anchors, in order on both sequences: level 0 over the whole pair, then inside every gap
    still longer than `budget` letters the next level, or with `stay` the same one until a gap yields nothing.
    A gap with more than `max_anchors` anchors is left to the next level; `lookback` caps how far back each
    anchor's chain search goes (0: no cap)."""
    var out = List[Anchor]()
    _refine(partner, query, 0, 0, levels, 0, costs, budget, stay, max_anchors, lookback, out)
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
    mut out: List[Anchor],
) raises:
    if level >= len(levels) or max(len(partner), len(query)) <= budget:
        return
    var current = levels[level]
    var found = current.select(find_anchors(partner, query, current.min_len, current.counted_to()))
    var picked = List[Anchor]()
    if len(found) <= max_anchors:
        picked = without_overlaps(chain(found, len(partner), len(query), costs, lookback).anchors)
    var next = level if stay else level + 1
    if len(picked) == 0:
        _refine(
            partner, query, ref_offset, query_offset, levels, level + 1, costs, budget, stay, max_anchors, lookback, out
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
