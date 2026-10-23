"""An MSA by anchored progressive placement, refined on its sum of pairs, with a copy-level start for tandem arrays.

The distinct sequences join one at a time, longest first. Each new sequence is placed by the exact profile DP of
`refine.place`, in a band centred on anchors that every sequence already placed votes for: through its own seeds
against the new one, the first source to anchor a letter deciding it (sporos's MUMs of 51, 31 and 21 letters,
none inside a tandem repeat at the top level, then MUMs unique within each gap between those, which sit on copies
carrying something of their own; then SMEMs of 31; then the optimal pairwise alignment); a letter is anchored where
more than half the votes from at least two sequences name one column, and the anchors kept are the heaviest chain
rising in letter and column. The MSA is then refined (`refine`, with group moves).

Where every sequence holds one clear tandem array (`units.tandem_layout`), a second start is built: the arrays
aligned copy by copy, the flanks on either side placed as above, refined alike; the start whose refined objective is
lower is kept; so is a third start, placed and refined with every distinct sequence weighing one before a weighted
refinement, where it is lower (the objective: the pairs' differing columns plus `gap_open` per gap a pair opens, over every pair of
the input sequences, so a distinct sequence weighs as many as spell it, as vntr-eval's pairs and pgrealign's alleles
do). Where the sequences
hold a tandem repeat, both refinements know its period (`refine`'s copy-aware moves).
"""

from karan import Aligner, Costs, Mode, Ties

from .backbone import Level, MUMS, SMEMS, backbone
from .chain import ChainCosts
from .graph import DELETION, INSERTION, cigar_runs
from .refine import objective, place, refine
from .units import Layout, tandem_layout, tandem_repeat

comptime GAP = UInt8(ord("-"))
comptime TRIM = 10
"""Letters an anchor gives back at each end before voting."""
comptime SHARE = 0.5


struct _Source(Copyable, Movable):
    """One seed source of the placement's cascade, or (`optimal`) the optimal pairwise alignment."""

    var optimal: Bool
    var levels: List[Level]
    var repeat_k: Int
    var repeat_distance: Int
    var mask_depth: Int
    var one_indel: Int

    def __init__(
        out self,
        optimal: Bool,
        var levels: List[Level],
        repeat_k: Int,
        repeat_distance: Int,
        depth: Int,
        one_indel: Int = -1,
    ):
        self.optimal = optimal
        self.levels = levels^
        self.repeat_k = repeat_k
        self.repeat_distance = repeat_distance
        self.mask_depth = depth
        self.one_indel = one_indel


def _cascade(one_indel: Int) -> List[_Source]:
    """The sources, best first. With `one_indel` >= 0, the anchors inside tandem repeats must leave their gap's
    length difference to one side (`backbone`'s `one_indel`), and the SMEMs are masked like the MUMs. With -2, no
    seeds: the optimal pairwise alignments alone vote."""
    var out = List[_Source]()
    if one_indel == -2:
        out.append(_Source(True, List[Level](), 0, 0, Int.MAX))
        return out^
    var mums = List[Level]()
    for k in [51, 31, 21]:
        mums.append(Level(MUMS, k))
    out.append(_Source(False, mums^, 12, 500, 1, one_indel))
    var smems = List[Level]()
    smems.append(Level(SMEMS, 31))
    if one_indel >= 0:
        out.append(_Source(False, smems^, 12, 500, 1, one_indel))
    else:
        out.append(_Source(False, smems^, 0, 0, Int.MAX))
    out.append(_Source(True, List[Level](), 0, 0, Int.MAX))
    return out^


def _partner_columns(
    source: _Source,
    partner: List[UInt8],
    q: List[UInt8],
    columns: List[Int],
    mut aligner: Aligner,
    two_piece: Costs,
    mut vote_c: List[List[Int]],
) raises:
    """The columns one placed sequence (its letters' `columns`) votes for the new sequence `q`'s letters."""
    if source.optimal:
        var found = aligner.align(Span(partner), Span(q), two_piece, Mode.GLOBAL, ties=Ties.LEFT)
        var runs = cigar_runs(found.cigar)
        var i = 0
        var j = 0
        for x in range(len(runs.ops)):
            for _ in range(runs.counts[x]):
                var op = runs.ops[x]
                if op == DELETION:
                    i += 1
                elif op == INSERTION:
                    j += 1
                else:
                    vote_c[j].append(columns[i])
                    i += 1
                    j += 1
        return
    var anchors = backbone(
        Span(partner),
        Span(q),
        source.levels,
        ChainCosts(substitution=0.1),
        0,
        True,
        repeat_k=source.repeat_k,
        repeat_distance=source.repeat_distance,
        mask_depth=source.mask_depth,
        one_indel=source.one_indel,
    )
    for a in anchors:
        var lo = min(TRIM, a.length // 2)
        for t in range(lo, a.length - lo):
            vote_c[a.query_pos + t].append(columns[a.ref_pos + t])


def _leading(
    source: _Source,
    q: List[UInt8],
    placed: List[List[UInt8]],
    column_of: List[List[Int]],
    mut aligner: Aligner,
    two_piece: Costs,
    alone: Bool,
) raises -> List[Int]:
    """Each letter of `q`'s column by `source`'s votes from every placed sequence (-1: none): more than SHARE of its
    votes, from at least two sequences (`alone`: a single placed sequence's vote suffices)."""
    var n = len(q)
    var vote_c = List[List[Int]]()
    for _ in range(n):
        vote_c.append(List[Int]())
    for r in range(len(placed)):
        if len(placed[r]) > 0:
            _partner_columns(source, placed[r], q, column_of[r], aligner, two_piece, vote_c)
    var leading = List[Int](length=n, fill=-1)
    for j in range(n):
        var total = len(vote_c[j])
        if total < (1 if alone else 2):
            continue
        sort(vote_c[j])
        var best_c = -1
        var best_w = 0
        var run = 0
        for x in range(total):
            run = run + 1 if x > 0 and vote_c[j][x] == vote_c[j][x - 1] else 1
            if run > best_w:
                best_w = run
                best_c = vote_c[j][x]
        if Float64(best_w) > SHARE * Float64(total):
            leading[j] = best_c
    return leading^


def _heaviest_chain(leading: List[Int], weight: List[Int], columns: Int) -> List[Int]:
    """The letters whose columns form the heaviest chain rising in letter and column (a Fenwick tree of the best
    chain ending before each column)."""
    var n = len(leading)
    var tree_w = List[Int](length=columns + 2, fill=0)
    var tree_i = List[Int](length=columns + 2, fill=-1)
    var best_w = List[Int](length=n, fill=0)
    var back = List[Int](length=n, fill=-1)
    var top = -1
    var top_w = 0
    for j in range(n):
        var c = leading[j]
        if c < 0:
            continue
        var w = 0
        var i = -1
        var k = c
        while k > 0:
            if tree_w[k] > w:
                w = tree_w[k]
                i = tree_i[k]
            k -= k & (-k)
        best_w[j] = w + weight[j]
        back[j] = i
        if best_w[j] > top_w:
            top_w = best_w[j]
            top = j
        k = c + 1
        while k <= columns + 1:
            if best_w[j] > tree_w[k]:
                tree_w[k] = best_w[j]
                tree_i[k] = j
            k += k & (-k)
    var chain = List[Int]()
    var t = top
    while t >= 0:
        chain.append(t)
        t = back[t]
    chain.reverse()
    return chain^


def _centre_of(points_c: List[Int], points_j: List[Int], columns: Int, n: Int) -> List[Int]:
    """Letters before each column, interpolated through (column, letters before it) points (rising)."""
    var pc = List[Int]()
    var pj = List[Int]()
    pc.append(0)
    pj.append(0)
    for t in range(len(points_c)):
        if points_c[t] > pc[len(pc) - 1] and points_j[t] >= pj[len(pj) - 1]:
            pc.append(points_c[t])
            pj.append(points_j[t])
            pc.append(points_c[t] + 1)
            pj.append(points_j[t] + 1)
    if pc[len(pc) - 1] < columns:
        pc.append(columns)
        pj.append(n)
    var centre = List[Int](length=columns + 1, fill=0)
    var p = 0
    for c in range(columns + 1):
        while p + 1 < len(pc) - 1 and pc[p + 1] <= c:
            p += 1
        var c0 = pc[p]
        var c1 = pc[min(p + 1, len(pc) - 1)]
        var j0 = pj[p]
        var j1 = pj[min(p + 1, len(pj) - 1)]
        centre[c] = min(n, max(0, j0 if c1 <= c0 else j0 + (c - c0) * (j1 - j0) // (c1 - c0)))
    return centre^


def _placed(
    seqs: List[List[UInt8]],
    band: Int,
    max_band: Int,
    gap_open: Int,
    weights: List[Int],
    one_indel: Int,
    vote_costs: Costs,
    mismatch: Int,
    gap: Int,
    hard: Int,
    nearest: Bool,
) raises -> List[List[UInt8]]:
    """The sequences (distinct or not, any may be empty) placed one at a time, longest first (with `nearest`: the
    longest, then each time the sequence nearest by edit distance to one already placed, the most confidently
    placed), each placement counting the placed rows' `weights`: rows in `seqs`' order. The letters the first `hard` sources anchor (and the heaviest
    chain keeps) are held to their columns, the rest placed by the DP between them; with 0, every anchor only centres
    the band."""
    var n = len(seqs)
    var keys = List[UInt64]()
    for k in range(n):
        keys.append((UInt64((1 << 32) - len(seqs[k])) << 32) | UInt64(k))
    sort(keys)
    var order = List[Int]()
    for key in keys:
        order.append(Int(key & 0xFFFFFFFF))
    if nearest and n > 2:
        # Prim's order from the longest: each next sequence the one nearest any placed one.
        var dist_aligner = Aligner()
        var near = List[Int](length=n, fill=Int.MAX)
        var done = List[Bool](length=n, fill=False)
        var prim = List[Int]()
        var at = order[0]
        for _ in range(n):
            done[at] = True
            prim.append(at)
            var next = -1
            for k in range(n):
                if done[k]:
                    continue
                var d = dist_aligner.distance(Span(seqs[at]), Span(seqs[k]))
                if d < near[k]:
                    near[k] = d
                if next < 0 or near[k] < near[next]:
                    next = k
            if next < 0:
                break
            at = next
        order = prim^
    var cascade = _cascade(one_indel)
    var aligner = Aligner()
    var two_piece = vote_costs
    var rows = List[List[UInt8]]()
    var placed = List[List[UInt8]]()
    var placed_weights = List[Int]()
    for step in range(n):
        var q = seqs[order[step]].copy()
        placed_weights.append(weights[order[step]])
        if step == 0:
            rows.append(q.copy())
            placed.append(q^)
            continue
        var columns = len(rows[0])
        var m = len(q)
        if m == 0:
            rows.append(List[UInt8](length=columns, fill=GAP))
            placed.append(q^)
            continue
        # The new row first takes columns of its own after the others' (every letter inserted).
        for r in range(len(rows)):
            for _ in range(m):
                rows[r].append(GAP)
        var mine = List[UInt8](length=columns, fill=GAP)
        mine.extend(q.copy())
        rows.append(mine^)
        var column_of = List[List[Int]]()
        for r in range(len(placed)):
            var cols = List[Int]()
            for c in range(len(rows[r])):
                if rows[r][c] != GAP:
                    cols.append(c)
            column_of.append(cols^)
        var leading = List[Int](length=m, fill=-1)
        var weight = List[Int](length=m, fill=0)
        var top_weight = 1
        for _ in range(len(cascade) - 1):
            top_weight *= 1000
        if step == 1:
            # A single placed sequence votes alone, through the first source's seeds.
            leading = _leading(cascade[0], q, placed, column_of, aligner, two_piece, True)
            for j in range(m):
                weight[j] = top_weight
        else:
            var rank_weight = 1
            for _ in range(len(cascade)):
                rank_weight *= 1000
            for source in cascade:
                rank_weight //= 1000
                var got = _leading(source, q, placed, column_of, aligner, two_piece, False)
                for j in range(m):
                    if leading[j] < 0 and got[j] >= 0:
                        leading[j] = got[j]
                        weight[j] = rank_weight
        placed.append(q^)
        var total_columns = len(rows[0])
        var chain_j = _heaviest_chain(leading, weight, total_columns)
        var pc = List[Int]()
        var pj = List[Int]()
        for j in chain_j:
            pc.append(leading[j])
            pj.append(j)
        var centre = _centre_of(pc, pj, total_columns, m)
        if len(chain_j) == 0:
            # No anchor: the diagonal through the others' columns.
            for c in range(total_columns + 1):
                centre[c] = min(m, c * m // max(1, columns))
        # Held letters: those the first `hard` sources anchored, on the chain.
        var forced = List[Int]()
        if hard > 0 and len(chain_j) > 0:
            var threshold = top_weight
            for _ in range(hard - 1):
                threshold //= 1000
            forced = List[Int](length=m, fill=-1)
            for j in chain_j:
                if weight[j] >= threshold:
                    forced[j] = leading[j]
        _ = place(rows, len(rows) - 1, centre, forced, band, max_band, gap_open, placed_weights, mismatch, gap)
    var out = List[List[UInt8]]()
    for _ in range(n):
        out.append(List[UInt8]())
    for t in range(n):
        out[order[t]] = rows[t].copy()
    return out^


def _letters(seq: List[UInt8], a: Int, b: Int) -> List[UInt8]:
    var out = List[UInt8](capacity=b - a)
    for k in range(a, b):
        out.append(seq[k])
    return out^


def progressive(
    seqs: List[List[UInt8]],
    passes: Int = 20,
    band: Int = 64,
    max_band: Int = 4096,
    gap_open: Int = 1,
    group_band: Int = 256,
    units: Bool = True,
    even_start: Bool = True,
    one_indel: Int = -1,
    unit_votes: Bool = False,
    mismatch: Int = 1,
    gap: Int = 1,
    hard: Int = 0,
    nearest: Bool = False,
) raises -> List[List[UInt8]]:
    """The sequences' MSA (module docstring): one row per sequence, in order, identical sequences sharing their row.
    `passes`, `band`, `max_band`, `gap_open` and `group_band` are the refinement's (see `refine`); `units` builds the
    copy-level start too where the sequences hold a tandem array, and gives the refinement the repeat's period;
    `even_start`, where some sequence repeats, one more start placed and refined with every distinct sequence
    weighing one, then refined weighted (a different path to a local optimum, which is sometimes lower);
    `one_indel`, the seeds' rule inside tandem repeats (`_cascade`); `unit_votes`, the optimal pairwise alignments that
    vote for placement and centre the refinement's bands at unit costs (vntr-eval's), not two-piece ones; `mismatch`
    and `gap`, the prices of a column of two different letters and of a letter against a gap (`refine`); `hard`, how
    many of the best seed sources' anchors hold their letters in placement (`_placed`); `nearest`, sequences placed
    nearest first, not longest first (`_placed`)."""
    var n = len(seqs)
    var distinct = List[List[UInt8]]()
    var weight = List[Int]()
    var of = List[Int](length=n, fill=-1)
    for r in range(n):
        for d in range(len(distinct)):
            if distinct[d] == seqs[r]:
                of[r] = d
                weight[d] += 1
                break
        if of[r] < 0:
            of[r] = len(distinct)
            distinct.append(seqs[r].copy())
            weight.append(1)
    var two_piece = Costs.edit() if unit_votes else Costs.two_piece(4, 4, 2, 24, 1)
    # The tandem repeat's copies, where the sequences hold one: the copy-level start, and the period the
    # refinement's copy-aware moves use (from the layout, else the longest sequence's array).
    var layout = tandem_layout(distinct, weights=weight) if units and len(distinct) > 2 else Layout()
    var period = layout.period
    if units and period == 0 and len(distinct) > 1:
        var longest = 0
        for k in range(1, len(distinct)):
            if len(distinct[k]) > len(distinct[longest]):
                longest = k
        var top = tandem_repeat(distinct[longest])
        if top.identity >= 0.7 and (top.end - top.start) * 10 >= len(distinct[longest]) * 4:
            period = top.period
    var rows = _placed(distinct, band, max_band, gap_open, weight, one_indel, two_piece, mismatch, gap, hard, nearest)
    if passes > 0 and len(distinct) > 1:
        rows = refine(
            rows,
            passes,
            band,
            max_band,
            True,
            two_piece,
            gap_open,
            False,
            group_band,
            period,
            weight,
            1000,
            mismatch,
            gap,
        )
    var repeated = len(distinct) < n
    if even_start and repeated and passes > 0 and len(distinct) > 1:
        var even = List[Int](length=len(distinct), fill=1)
        var other = _placed(
            distinct, band, max_band, gap_open, even, one_indel, two_piece, mismatch, gap, hard, nearest
        )
        other = refine(
            other,
            passes,
            band,
            max_band,
            True,
            two_piece,
            gap_open,
            False,
            group_band,
            period,
            even,
            1000,
            mismatch,
            gap,
        )
        other = refine(
            other,
            passes,
            band,
            max_band,
            True,
            two_piece,
            gap_open,
            False,
            group_band,
            period,
            weight,
            1000,
            mismatch,
            gap,
        )
        if objective(other, gap_open, weight, mismatch, gap) < objective(rows, gap_open, weight, mismatch, gap):
            rows = other^
    if units and len(distinct) > 2:
        if len(layout.rows) == len(distinct):
            var left = List[List[UInt8]]()
            var right = List[List[UInt8]]()
            for k in range(len(distinct)):
                left.append(_letters(distinct[k], 0, layout.starts[k]))
                right.append(_letters(distinct[k], layout.ends[k], len(distinct[k])))
            var left_rows = _placed(
                left, band, max_band, gap_open, weight, one_indel, two_piece, mismatch, gap, hard, nearest
            )
            var right_rows = _placed(
                right, band, max_band, gap_open, weight, one_indel, two_piece, mismatch, gap, hard, nearest
            )
            var other = List[List[UInt8]]()
            for k in range(len(distinct)):
                var row = left_rows[k].copy()
                row.extend(layout.rows[k].copy())
                row.extend(right_rows[k].copy())
                other.append(row^)
            if passes > 0:
                other = refine(
                    other,
                    passes,
                    band,
                    max_band,
                    True,
                    two_piece,
                    gap_open,
                    False,
                    group_band,
                    period,
                    weight,
                    1000,
                    mismatch,
                    gap,
                )
            if objective(other, gap_open, weight, mismatch, gap) < objective(rows, gap_open, weight, mismatch, gap):
                rows = other^
    var out = List[List[UInt8]]()
    for r in range(n):
        out.append(rows[of[r]].copy())
    return out^
