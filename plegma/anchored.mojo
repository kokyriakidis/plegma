"""Anchored-read quality checks, foreign-read detection, and one-directional reach/trimming for `poa.mojo`.

Separated from `poa.mojo` so end-checking (`ends_checked`), paralog/other-copy detection (`find_foreign`), and
one-sided reach/trimming (`read_reach`, `trim_consensus`) can be tested and parallelized independently.
"""

from std.atomic import Atomic

from karan import Aligner, Alignment, Costs, Mode
from karan.common import spread, thread_count

from .bands import past_the_end
from .graph import Graph
from .kernel import BOTH, LEFT, RIGHT

comptime FOREIGN_SPANNING = 3
"""The fewest spanning reads whose consensus a read anchored at one end is checked against (`find_foreign`)."""
comptime FOREIGN_RATIO = 1.5
comptime FOREIGN_SLACK = 3
"""A read anchored at one end is another copy's when its edit distance to the spanning reads' consensus passes
`FOREIGN_RATIO` times the worst spanning read's rate over the letters it covers, plus `FOREIGN_SLACK`."""

comptime END_LETTERS = 100
"""A spanning read's end checked against the reads anchored at that end (`ends_checked`): its first or last this
many letters."""
comptime END_EDITS = 35
"""Two reads' ends agree within this many edits over `END_LETTERS`."""
comptime _STRAY = -1
"""A read given as spanning whose ends both disagree with the reads anchored there (`ends_checked`)."""


def _ends_agree(mut aligner: Aligner, a: List[UInt8], b: List[UInt8], at_start: Bool) raises -> Bool:
    """Whether `a`'s and `b`'s first (`at_start`) or last `END_LETTERS` letters are within `END_EDITS` edits."""
    var x = Span(a)[0:END_LETTERS] if at_start else Span(a)[len(a) - END_LETTERS :]
    var y = Span(b)[0:END_LETTERS] if at_start else Span(b)[len(b) - END_LETTERS :]
    return aligner.align(x, y, Costs.edit(), Mode.GLOBAL).cost <= END_EDITS


def _check_one_spanning(
    mut aligner: Aligner,
    seqs: List[List[UInt8]],
    s: Int,
    spanning: List[Int],
    left: List[Int],
    right: List[Int],
) raises -> Int:
    """Determines the verified anchor (`BOTH`, `LEFT`, `RIGHT`, or `_STRAY`) for spanning read `s`."""
    var held = 0
    for side in range(2):
        var at_start = side == 0
        ref ones = left if at_start else right
        if len(ones) == 0:
            held += 1 if at_start else 2
            continue
        var shared = False
        for r in spanning:
            if r != s and _ends_agree(aligner, seqs[s], seqs[r], at_start):
                shared = True
                break
        if shared:
            held += 1 if at_start else 2
            continue
        var agree = 0
        var disagree = 0
        for r in ones:
            if r == s or len(seqs[r]) < END_LETTERS:
                continue
            if _ends_agree(aligner, seqs[s], seqs[r], at_start):
                agree += 1
            else:
                disagree += 1
        if disagree <= agree:
            held += 1 if at_start else 2
    if held == 0:
        return _STRAY
    if held == 1:
        return LEFT
    if held == 2:
        return RIGHT
    return BOTH


def ends_checked(
    mut aligner: Aligner, seqs: List[List[UInt8]], anchors: List[Int], threads: Int = 1
) raises -> List[Int]:
    """`anchors`, each spanning read's ends checked against the reads anchored at that end."""
    var out = anchors.copy()
    var left = List[Int]()
    var right = List[Int]()
    var spanning = List[Int]()
    for r in range(len(seqs)):
        if anchors[r] == LEFT:
            left.append(r)
        elif anchors[r] == RIGHT:
            right.append(r)
        elif len(seqs[r]) >= END_LETTERS:
            spanning.append(r)
    if len(spanning) == 0 or len(left) + len(right) == 0:
        return out^
    var workers = thread_count(threads, len(spanning))
    if workers > 1 and len(spanning) >= 2:
        var results = List[Int](length=len(spanning), fill=BOTH)
        var failed = List[Bool](length=len(spanning), fill=False)
        var results_p = results.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def check_worker(
            slot: Int,
        ) {mut taken, imm spanning, imm seqs, imm left, imm right, imm results_p, imm failed_p}:
            var local_aligner = Aligner(warm=False)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= len(spanning):
                    return
                try:
                    results_p[unsafe_offset=idx] = _check_one_spanning(
                        local_aligner, seqs, spanning[idx], spanning, left, right
                    )
                except:
                    failed_p[unsafe_offset=idx] = True

        spread(check_worker, workers, workers)
        for idx in range(len(spanning)):
            var s = spanning[idx]
            var anchor = _check_one_spanning(aligner, seqs, s, spanning, left, right) if failed[idx] else results[idx]
            if anchor != BOTH:
                out[s] = anchor
        return out^
    for s in spanning:
        var anchor = _check_one_spanning(aligner, seqs, s, spanning, left, right)
        if anchor != BOTH:
            out[s] = anchor
    return out^


def as_built(anchor_of: List[Int]) -> List[Int]:
    """Each read's anchor as `poa` built with it, a stray as `BOTH` (it is in `Poa.foreign`, and polishes nothing)."""
    var out = anchor_of.copy()
    for r in range(len(out)):
        if out[r] == _STRAY:
            out[r] = BOTH
    return out^


def _one_sided_align(mut aligner: Aligner, consensus: List[UInt8], seq: List[UInt8], left: Bool) raises -> Alignment:
    """`seq`, anchored at the consensus's start (`left`) or end, against it: ending inside it or past it, the cheaper.
    """
    var found = aligner.align(Span(consensus), Span(seq), Costs.edit(), Mode.PREFIX if left else Mode.SUFFIX)
    if len(consensus) - len(seq) < found.cost:
        var other = past_the_end(aligner, consensus, seq, left)
        if other.cost < found.cost:
            found = other^
    return found^


def _far(mut aligner: Aligner, core: List[UInt8], seq: List[UInt8], left: Bool, worst: Float64) raises -> Bool:
    """Whether `seq`, anchored at `core`'s start (`left`) or end, is further from it than `find_foreign` allows."""
    var found = _one_sided_align(aligner, core, seq, left)
    var covered = found.reference_end - found.reference_start
    return Float64(found.cost) > FOREIGN_RATIO * worst * Float64(covered) + Float64(FOREIGN_SLACK)


def find_foreign(
    mut aligner: Aligner,
    core: List[UInt8],
    seqs: List[List[UInt8]],
    spanning: List[Int],
    left: List[Int],
    right: List[Int],
    threads: Int = 1,
) raises -> List[Int]:
    """The reads of `left` and `right` another copy's: further from `core`, the spanning reads' consensus, than
    `FOREIGN_RATIO` times the worst spanning read's rate over the letters each covers, plus `FOREIGN_SLACK`."""
    var worst = 0.0
    var span_workers = thread_count(threads, len(spanning))
    if span_workers > 1 and len(spanning) >= 2:
        var rates = List[Float64](length=len(spanning), fill=0.0)
        var failed = List[Bool](length=len(spanning), fill=False)
        var rates_p = rates.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)
        var denom = Float64(max(1, len(core)))

        def span_worker(slot: Int) {mut taken, imm spanning, imm core, imm seqs, imm denom, imm rates_p, imm failed_p}:
            var local_aligner = Aligner(warm=False)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= len(spanning):
                    return
                try:
                    var found = local_aligner.align(Span(core), Span(seqs[spanning[idx]]), Costs.edit(), Mode.GLOBAL)
                    rates_p[unsafe_offset=idx] = Float64(found.cost) / denom
                except:
                    failed_p[unsafe_offset=idx] = True

        spread(span_worker, span_workers, span_workers)
        for idx in range(len(spanning)):
            if failed[idx]:
                var found = aligner.align(Span(core), Span(seqs[spanning[idx]]), Costs.edit(), Mode.GLOBAL)
                worst = max(worst, Float64(found.cost) / denom)
            else:
                worst = max(worst, rates[idx])
    else:
        for r in spanning:
            var found = aligner.align(Span(core), Span(seqs[r]), Costs.edit(), Mode.GLOBAL)
            worst = max(worst, Float64(found.cost) / Float64(max(1, len(core))))

    var total_ones = len(left) + len(right)
    var out = List[Int]()
    var one_workers = thread_count(threads, total_ones)
    if one_workers > 1 and total_ones >= 2:
        var is_far = List[Bool](length=total_ones, fill=False)
        var failed = List[Bool](length=total_ones, fill=False)
        var far_p = is_far.unsafe_ptr()
        var failed_p = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)
        var n_left = len(left)

        def ones_worker(
            slot: Int,
        ) {
            mut taken,
            imm total_ones,
            imm n_left,
            imm left,
            imm right,
            imm core,
            imm seqs,
            imm worst,
            imm far_p,
            imm failed_p,
        }:
            var local_aligner = Aligner(warm=False)
            while True:
                var idx = Int(taken.fetch_add(1))
                if idx >= total_ones:
                    return
                var is_left = idx < n_left
                var r = left[idx] if is_left else right[idx - n_left]
                try:
                    far_p[unsafe_offset=idx] = _far(local_aligner, core, seqs[r], is_left, worst)
                except:
                    failed_p[unsafe_offset=idx] = True

        spread(ones_worker, one_workers, one_workers)
        for idx in range(total_ones):
            var is_left = idx < n_left
            var r = left[idx] if is_left else right[idx - n_left]
            var far = _far(aligner, core, seqs[r], is_left, worst) if failed[idx] else is_far[idx]
            if far:
                out.append(r)
    else:
        for r in left:
            if _far(aligner, core, seqs[r], True, worst):
                out.append(r)
        for r in right:
            if _far(aligner, core, seqs[r], False, worst):
                out.append(r)
    sort(out)
    return out^


def without_reads(reads: List[Int], drop: List[Int]) -> List[Int]:
    """`reads` but those in `drop` (sorted)."""
    var out = List[Int]()
    for r in reads:
        var lo = 0
        var hi = len(drop)
        while lo < hi:
            var mid = (lo + hi) // 2
            if drop[mid] < r:
                lo = mid + 1
            else:
                hi = mid
        if lo == len(drop) or drop[lo] != r:
            out.append(r)
    return out^


def read_reach(graph: Graph, consensus: List[Int], paths: List[Int], from_left: Bool) -> List[Int]:
    """How far along `consensus` (its nodes) each of `paths` reaches from the anchor it starts at."""
    var at_column = List[Int](length=graph.column_count(), fill=-1)
    for k in range(len(consensus)):
        at_column[Int(graph.columns[consensus[k]])] = k
    var out = List[Int](capacity=len(paths))
    for p in paths:
        var furthest = 0
        for node in graph.paths[p]:
            var k = at_column[Int(graph.columns[Int(node)])]
            if k >= 0:
                furthest = max(furthest, k + 1 if from_left else len(consensus) - k)
        out.append(furthest)
    return out^


def trim_consensus(mut consensus: List[UInt8], reach: List[Int], least: Int, from_left: Bool):
    """A one-directional consensus cut where fewer than `least` reads go on."""
    if least <= 0 or len(reach) == 0:
        return
    var reaches = reach.copy()
    sort(reaches)
    var cut = reaches[len(reaches) - min(least, len(reaches))]
    if cut >= len(consensus):
        return
    if from_left:
        consensus.resize(cut, 0)
    else:
        var kept = List[UInt8](capacity=cut)
        for i in range(len(consensus) - cut, len(consensus)):
            kept.append(consensus[i])
        consensus = kept^
