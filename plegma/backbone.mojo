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
    var outside_repeats: Bool
    """Its anchors kept out of the partner's tandem repeats (`backbone`'s mask) at every depth: the flanks."""

    def __init__(out self, kind: Int, min_len: Int, cap: Int = 1, outside_repeats: Bool = False):
        self.kind = kind
        self.min_len = min_len
        self.cap = cap
        self.outside_repeats = outside_repeats

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
    one_indel: Int = -1,
) raises -> List[Anchor]:
    """Non-overlapping anchors, in order on both sequences: level 0 over the whole pair, then inside every gap
    still longer than `budget` letters the next level, or with `stay` the same one until a gap yields nothing.
    A gap with more than `max_anchors` anchors is left to the next level; `lookback` caps how far back each
    anchor's chain search goes (0: no cap). With `repeat_k`, the anchors inside tandem repeats of the partner
    (`RepeatMask` of `repeat_k`-mers within `repeat_distance`) are dropped before chaining, at the first
    `mask_depth` depths of the recursion (every depth by default): with 1, the whole pair's anchors keep out of
    the repeats and so pin the flanks around them, and inside a gap between anchors, a repeat array, a match is
    unique within the gap's two pieces, as a copy that carries something of its own is. With `one_indel` >= 0,
    past the mask's depths an anchor inside the repeats is kept only if its diagonal lies within that many letters
    of one of its gap's ends: the gap's length difference (the copies one piece has and the other lacks) all before
    it or all after it. An optimal alignment of two arrays takes that difference as one indel, and an anchor
    splitting it pairs copies no optimal alignment pairs (on vntr-eval's regions, three anchored letters in four
    of those, against one in fifteen of the others)."""
    var out = List[Anchor]()
    var repeats = RepeatMask(partner, repeat_k, repeat_distance)
    _refine(
        partner,
        query,
        0,
        0,
        levels,
        0,
        costs,
        budget,
        stay,
        max_anchors,
        lookback,
        repeats,
        mask_depth,
        one_indel,
        0,
        out,
    )
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
    one_indel: Int,
    depth: Int,
    mut out: List[Anchor],
) raises:
    if level >= len(levels) or max(len(partner), len(query)) <= budget:
        return
    var current = levels[level]
    var found = current.select(find_anchors(partner, query, current.min_len, current.counted_to()))
    if repeats.k > 0 and (depth < mask_depth or current.outside_repeats):
        var outside = List[Anchor]()
        for a in found:
            if not repeats.covers(ref_offset + a.ref_pos, ref_offset + a.ref_end()):
                outside.append(a)
        found = outside^
    elif repeats.k > 0 and one_indel >= 0:
        var kept = List[Anchor]()
        for a in found:
            var start = a.query_pos - a.ref_pos
            var end = (len(query) - a.query_end()) - (len(partner) - a.ref_end())
            if (
                not repeats.covers(ref_offset + a.ref_pos, ref_offset + a.ref_end())
                or abs(start) <= one_indel
                or abs(end) <= one_indel
            ):
                kept.append(a)
        found = kept^
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
            one_indel,
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
            one_indel,
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
        one_indel,
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


def split_anchors(
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    min_len: Int = 51,
    len_d: Float64 = 0.5,
    max_depth: Int = 3,
    one_indel: Int = 2,
    repeat_k: Int = 12,
    repeat_distance: Int = 500,
    ambiguity: Float64 = 0.1,
) raises -> List[Anchor]:
    """Non-overlapping anchors in order on both sequences, placed where they are sure. First the flanks: MUMs of 51,
    31 and 21 letters outside the partner's tandem repeats, chained (`backbone`), which pin where the alignment
    starts and ends. Then, in every window they leave, one anchor at a time: the longest match unique in both of the
    window's pieces (its frequencies counted in the window alone), of `min_len` letters or more, and of `len_d` times
    the window's length difference D or more (the larger the gap a window must place, the more an anchor must hold
    to decide where it goes), that leaves D to one side of it (its diagonal within `one_indel` letters of one end:
    one indel); then the same in the two windows it leaves, `max_depth` deep at most. An anchor with a copy
    (`_ambiguous`: its letters on another diagonal of its window with at most `ambiguity` of them mismatched) is
    dropped from the flanks and passed over in the windows: inside tandem arrays a long exact match unique in both
    sequences is often a copy the optimum pairs elsewhere, as in repeats the mask misses (copies over
    `repeat_distance` apart, or diverged) and in windows whose length difference is no single indel. On vntr-eval's
    regions (426 pairs, every anchor checked against every optimal alignment, `anchors --brute`) the anchors force
    no excess over the optima under gap-affine prices (vntr-eval's 4, 6 + 2k and plegma's 2, 1 + 2k), against
    0.04-0.06% without the copy test and 2-6% for chained MUMs, SMEMs or MEMs; under unit cost 0.01%, ties where a
    long deletion spreads into many small ones."""
    var out = List[Anchor]()
    var flank_levels = List[Level]()
    for k in [51, 31, 21]:
        flank_levels.append(Level(MUMS, k, 1, True))
    var flanks = backbone(
        partner,
        query,
        flank_levels,
        ChainCosts(substitution=0.1),
        0,
        True,
        repeat_k=repeat_k,
        repeat_distance=repeat_distance,
        mask_depth=0,
    )
    if ambiguity > 0.0:
        var kept = List[Anchor]()
        for f in flanks:
            if not _ambiguous(partner, query, 0, len(partner), 0, len(query), f, ambiguity):
                kept.append(f)
        flanks = kept^
    var ref_at = 0
    var query_at = 0
    for f in flanks:
        _split(
            partner,
            query,
            ref_at,
            f.ref_pos,
            query_at,
            f.query_pos,
            min_len,
            len_d,
            max_depth,
            one_indel,
            1,
            out,
            ambiguity,
        )
        out.append(f)
        ref_at = f.ref_end()
        query_at = f.query_end()
    _split(
        partner,
        query,
        ref_at,
        len(partner),
        query_at,
        len(query),
        min_len,
        len_d,
        max_depth,
        one_indel,
        1,
        out,
        ambiguity,
    )
    return out^


def _ambiguous(
    partner: ImmSpan[UInt8, _], query: ImmSpan[UInt8, _], r0: Int, r1: Int, q0: Int, q1: Int, x: Anchor, share: Float64
) -> Bool:
    """Whether anchor `x` has a copy in its window (`r0 ..< r1` of the partner, `q0 ..< q1` of the query): its
    partner letters against the query window, or its query letters against the partner window, on another diagonal
    with at most `share` of them mismatched (no indels)."""
    var limit = Int(share * Float64(x.length))
    for s in range(q0, q1 - x.length + 1):
        if s == x.query_pos:
            continue
        var mismatches = 0
        for k in range(x.length):
            if partner[x.ref_pos + k] != query[s + k]:
                mismatches += 1
                if mismatches > limit:
                    break
        if mismatches <= limit:
            return True
    for s in range(r0, r1 - x.length + 1):
        if s == x.ref_pos:
            continue
        var mismatches = 0
        for k in range(x.length):
            if query[x.query_pos + k] != partner[s + k]:
                mismatches += 1
                if mismatches > limit:
                    break
        if mismatches <= limit:
            return True
    return False


def _split(
    partner: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    r0: Int,
    r1: Int,
    q0: Int,
    q1: Int,
    min_len: Int,
    len_d: Float64,
    max_depth: Int,
    one_indel: Int,
    depth: Int,
    mut out: List[Anchor],
    ambiguity: Float64 = 0.0,
) raises:
    if r1 - r0 < min_len or q1 - q0 < min_len or depth > max_depth:
        return
    var d = (q1 - q0) - (r1 - r0)
    var need = max(min_len, Int(len_d * Float64(abs(d))))
    var found = mums(find_anchors(partner[r0:r1], query[q0:q1], need, 2))
    # The candidates leaving D to one side, longest first; with `ambiguity`, the longest without a copy.
    var fits = List[Int]()
    for k in range(len(found)):
        var start = found[k].query_pos - found[k].ref_pos
        if start < min(0, d) - one_indel or start > max(0, d) + one_indel:
            continue
        if abs(start) > one_indel and abs(d - start) > one_indel:
            continue
        fits.append(k)
    var best = -1
    while len(fits) > 0:
        var top = 0
        for u in range(1, len(fits)):
            if found[fits[u]].length > found[fits[top]].length:
                top = u
        var k = fits[top]
        if ambiguity <= 0.0 or not _ambiguous(partner, query, r0, r1, q0, q1, found[k].shifted(r0, q0), ambiguity):
            best = k
            break
        _ = fits.pop(top)
    if best < 0:
        return
    var x = found[best]
    _split(
        partner,
        query,
        r0,
        r0 + x.ref_pos,
        q0,
        q0 + x.query_pos,
        min_len,
        len_d,
        max_depth,
        one_indel,
        depth + 1,
        out,
        ambiguity,
    )
    out.append(x.shifted(r0, q0))
    _split(
        partner,
        query,
        r0 + x.ref_end(),
        r1,
        q0 + x.query_end(),
        q1,
        min_len,
        len_d,
        max_depth,
        one_indel,
        depth + 1,
        out,
        ambiguity,
    )
