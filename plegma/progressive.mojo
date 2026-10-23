"""An MSA by anchored progressive placement, refined on its sum of pairs, with a copy-level start for tandem arrays.

The distinct sequences join one at a time, longest first. Each new sequence is placed by the exact profile DP of
`refine.place`, in a band centred on anchors that every sequence already placed votes for: by default each placed
sequence's optimal unit-cost alignment to the new one names a column for each of its letters; a letter is anchored
where more than half the votes from at least two sequences name one column, and the anchors kept are the heaviest
chain rising in letter and column. Seeds can vote first (`one_indel` >= -1: sporos's MUMs of 51, 31 and 21 letters
outside tandem repeats, then inside each array, then SMEMs of 31, the alignment voting only where they do not), but
on vntr-eval's regions the anchors seeds give inside repeats pair letters no optimal alignment pairs (the anchor
audit, experiments/anchors: MUMs unique in a repeat window force 1.8% more than the optimum, two-piece alignments
12%, unit-cost alignments nothing), and placement from the alignments alone is as good or better and faster.
Anchors centre the band only: holding their letters to their columns locks in the placed rows' own errors.

The MSA is then refined (`refine`, with group moves) on the sum of pairs over every pair of the input sequences, a
distinct sequence weighing as many as spell it (as vntr-eval's pairs and pgrealign's alleles do), a column of two
letters costing `mismatch` (2), of a letter and a gap `gap` (2), and each gap a pair opens `gap_open` (1): half a
mismatch, where vntr-eval's unit and affine costs (0 and 1.5) trade off best (its medians beaten on both in every
stratum). Where the sequences hold a tandem repeat, the refinement knows its period (`refine`'s copy-aware moves).

Where every sequence holds one clear tandem array (`units.tandem_layout`), a second start is built: the arrays
aligned copy by copy, the flanks on either side placed as above, refined alike; and where a sequence repeats, a
third, placed and refined with every distinct sequence weighing one before a weighted refinement; the start whose
refined objective is lowest is kept.
"""

from std.os import getenv
from std.sys import stderr

from karan import Aligner, Costs, Mode, Ties

from .backbone import Level, MUMS, SMEMS, backbone, split_anchors
from .chain import ChainCosts
from .graph import DELETION, INSERTION, cigar_runs
from .refine import objective, place, refine, row_cost
from .units import Layout, tandem_layout, tandem_repeat

comptime TWO_PIECE_VOTES = 0
"""`progressive`'s votes at minimap2's two-piece costs (-O4,24 -E2,1, mismatch 4)."""
comptime UNIT_VOTES = 1
"""`progressive`'s votes at unit costs (vntr-eval's unit metric)."""
comptime PRICED_VOTES = 2
"""`progressive`'s votes gap-affine at its objective's prices."""
comptime AFFINE_VOTES = 3
"""`progressive`'s votes at vntr-eval's affine costs (mismatch 4, a gap of k letters 6 + 2k)."""
comptime GAP = UInt8(ord("-"))
comptime TRIM = 10
"""Letters an anchor gives back at each end before voting."""
comptime SHARE = 0.5
comptime SUBFAMILY_ROWS = 64
"""Up to this many distinct sequences, the subfamily start is built too (its tree costs their pairs' distances)."""
comptime VOTERS = 24
"""Past this many placed sequences, only the nearest this many (by edit distance) vote for a new one's letters: the
votes need a few dozen, and aligning to every placed sequence costs the square."""


struct _Source(Copyable, Movable):
    """One seed source of the placement's cascade, or (`optimal`) the optimal pairwise alignment."""

    var optimal: Bool
    var levels: List[Level]
    var repeat_k: Int
    var repeat_distance: Int
    var mask_depth: Int
    var one_indel: Int
    var split: Bool
    """The anchors `split_anchors` places, not a backbone's."""

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
        self.split = False


def _staged_levels() -> List[Level]:
    """The staged seed levels: MUMs of 51, 31 and 21 letters outside the repeats, then MUMs of 1,000, 600 and 400."""
    var levels = List[Level]()
    for k in [51, 31, 21]:
        levels.append(Level(MUMS, k, 1, True))
    for k in [1000, 600, 400]:
        levels.append(Level(MUMS, k))
    return levels^


def _cascade(one_indel: Int) -> List[_Source]:
    """The sources, best first. With `one_indel` >= 0, the anchors inside tandem repeats must leave their gap's
    length difference to one side (`backbone`'s `one_indel`), and the SMEMs are masked like the MUMs. With -2, no
    seeds: the optimal pairwise alignments alone vote; with -3, the longest SMEMs (LONG_SMEM letters or more) first,
    then the optimal alignments; with -4, the staged seeds (`_staged_levels`), then the optimal alignments; with -5,
    the split seeds (`backbone.split_anchors`), then the optimal alignments. A pair too long to align exactly
    (OPTIMAL_CELLS) votes through the split seeds in any case."""
    var out = List[_Source]()
    if one_indel == -2:
        out.append(_Source(True, List[Level](), 0, 0, Int.MAX))
        return out^
    if one_indel == -5:
        # The split seeds (`backbone.split_anchors`), then the optimal alignments.
        var seeds = _Source(False, List[Level](), 0, 0, Int.MAX)
        seeds.split = True
        out.append(seeds^)
        out.append(_Source(True, List[Level](), 0, 0, Int.MAX))
        return out^
    if one_indel == -4:
        # Staged seeds: the flanks first (MUMs outside the repeats, which pin where the alignment starts and ends),
        # then the longest unique matches (MUMs of 1,000, 600 and 400 letters) inside the gaps they leave; then the
        # optimal alignments (the anchor audit: a third the cost of MUMs unique in a window with the one-indel rule,
        # anchoring three and a half times the letters).
        out.append(_Source(False, _staged_levels(), 12, 500, 0))
        out.append(_Source(True, List[Level](), 0, 0, Int.MAX))
        return out^
    if one_indel == -3:
        # The longest SMEMs alone, then the optimal alignments: an exact match of LONG_SMEM letters reproduces a
        # long run of a repeat's copy variants in order, which a wrong copy rarely does (the anchor audit: SMEMs of
        # 600 or more force 0.1% over the optimum, against 3.8% for those of 31).
        var long_smems = List[Level]()
        long_smems.append(Level(SMEMS, LONG_SMEM))
        out.append(_Source(False, long_smems^, 0, 0, Int.MAX))
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
    if source.optimal and len(partner) * len(q) > OPTIMAL_CELLS:
        # Too long to align exactly: the split seeds vote instead.
        var seeds = _Source(False, List[Level](), 0, 0, Int.MAX)
        seeds.split = True
        _partner_columns(seeds, partner, q, columns, aligner, two_piece, vote_c)
        return
    if source.split:
        for a in split_anchors(Span(partner), Span(q)):
            for t in range(a.length):
                vote_c[a.query_pos + t].append(columns[a.ref_pos + t])
        return
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
    voters: List[Bool],
    mut support: List[Int],
    plurality: Bool = False,
) raises -> List[Int]:
    """Each letter of `q`'s column by `source`'s votes from the placed sequences `voters` marks (-1: none): more than
    SHARE of its votes (with `plurality`, the most), from at least two sequences (`alone`: a single placed sequence's
    vote suffices); `support[j]` the votes its column won."""
    var n = len(q)
    var vote_c = List[List[Int]]()
    for _ in range(n):
        vote_c.append(List[Int]())
    for r in range(len(placed)):
        if len(placed[r]) > 0 and voters[r]:
            _partner_columns(source, placed[r], q, column_of[r], aligner, two_piece, vote_c)
    var leading = List[Int](length=n, fill=-1)
    support = List[Int](length=n, fill=0)
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
        if plurality or Float64(best_w) > SHARE * Float64(total):
            leading[j] = best_c
            support[j] = best_w
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


comptime EDGE_SLACK = 8
comptime OPTIMAL_CELLS = 400_000_000
"""Past this many cells (a pair's lengths' product), a placed sequence votes through the staged seeds, not its
optimal alignment."""
comptime LONG_SMEM = 600
comptime CORRIDOR_CELLS = 150_000_000
"""Past this many cells, a pair's corridor is not computed (its forward costs are kept whole), and placement is
unconstrained."""


def _corridor(guide: List[UInt8], q: List[UInt8], slack: Int) -> Tuple[List[Int], List[Int]]:
    """The co-optimal corridor of `q` against `guide` at unit costs: for every count i of the guide's letters, the
    least and most of `q`'s letters an alignment at most `slack` over the optimum has consumed with them (a forward
    and a backward DP). Empty past CORRIDOR_CELLS."""
    var n = len(guide)
    var m = len(q)
    var lo = List[Int]()
    var hi = List[Int]()
    if (n + 1) * (m + 1) > CORRIDOR_CELLS:
        return (lo^, hi^)
    var w = m + 1
    var fwd = List[Int32](length=(n + 1) * w, fill=0)
    for j in range(w):
        fwd[j] = Int32(j)
    for i in range(1, n + 1):
        fwd[i * w] = Int32(i)
        var x = guide[i - 1]
        for j in range(1, w):
            var d = fwd[(i - 1) * w + j - 1] + (Int32(0) if x == q[j - 1] else Int32(1))
            fwd[i * w + j] = min(d, min(fwd[(i - 1) * w + j] + 1, fwd[i * w + j - 1] + 1))
    var best = fwd[n * w + m] + Int32(slack)
    lo = List[Int](length=n + 1, fill=0)
    hi = List[Int](length=n + 1, fill=m)
    var row = List[Int32](length=w, fill=0)
    var next_row = List[Int32](length=w, fill=0)
    for i in reversed(range(n + 1)):
        if i == n:
            for j in range(w):
                row[j] = Int32(m - j)
        else:
            row[m] = Int32(n - i)
            var x = guide[i]
            for j in reversed(range(m)):
                var d = next_row[j + 1] + (Int32(0) if x == q[j] else Int32(1))
                row[j] = min(d, min(next_row[j] + 1, row[j + 1] + 1))
        var first = -1
        var last = -1
        for j in range(w):
            if fwd[i * w + j] + row[j] <= best:
                if first < 0:
                    first = j
                last = j
        lo[i] = first if first >= 0 else 0
        hi[i] = last if last >= 0 else m
        next_row = row.copy()
    return (lo^, hi^)


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
    consistency: Bool,
    shuffle: Int = 0,
    exact_cells: Int = 0,
    interleave: Int = 0,
    corridor: Int = -1,
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
    if shuffle > 0 and n > 2:
        # A random order after the longest (a Fisher-Yates shuffle seeded by `shuffle`): a restart of the search.
        var state = UInt64(shuffle) * 6364136223846793005 + 1442695040888963407
        for t in reversed(range(2, n)):
            state = state * 6364136223846793005 + 1442695040888963407
            var u = 1 + Int((state >> 33) % UInt64(t))
            var keep = order[t]
            order[t] = order[u]
            order[u] = keep
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
    var auditing = getenv("PLEGMA_AUDIT") == "1"
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
        # The voters: every placed sequence, or past VOTERS the nearest that many (PLEGMA_VOTERS overrides it, for
        # experiments).
        var voter_cap = VOTERS
        var cap_text = getenv("PLEGMA_VOTERS")
        if cap_text != "":
            voter_cap = Int(cap_text)
        var voters = List[Bool](length=len(placed), fill=True)
        if len(placed) > voter_cap:
            var keys = List[UInt64]()
            for r in range(len(placed)):
                var d = aligner.distance(Span(placed[r]), Span(q)) if len(placed[r]) > 0 else Int(UInt32.MAX)
                keys.append((UInt64(d) << 32) | UInt64(r))
            sort(keys)
            for r in range(len(placed)):
                voters[r] = False
            for t in range(voter_cap):
                voters[Int(keys[t] & 0xFFFFFFFF)] = True
        var top_weight = 1
        for _ in range(len(cascade) - 1):
            top_weight *= 1000
        if step == 1:
            # A single placed sequence votes alone, through the first source's seeds.
            var support = List[Int]()
            leading = _leading(cascade[0], q, placed, column_of, aligner, two_piece, True, voters, support)
            for j in range(m):
                weight[j] = top_weight
        else:
            var rank_weight = 1
            for _ in range(len(cascade)):
                rank_weight *= 1000
            for source in cascade:
                rank_weight //= 1000
                var support = List[Int]()
                var got = _leading(
                    source, q, placed, column_of, aligner, two_piece, False, voters, support, consistency
                )
                for j in range(m):
                    if leading[j] < 0 and got[j] >= 0:
                        leading[j] = got[j]
                        # With `consistency`, an anchor weighs the votes its column won.
                        weight[j] = rank_weight * (support[j] if consistency else 1)
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
        if auditing:
            # The addition audit (PLEGMA_AUDIT): this placement's row cost, the best row possible against the placed
            # rows (the same DP unbanded), and the floor (its optimal pairwise costs to them at the same prices).
            var exact_rows = rows.copy()
            var unbanded = max(len(rows[0]), m) + 1
            _ = place(
                exact_rows,
                len(rows) - 1,
                List[Int](),
                List[Int](),
                unbanded,
                unbanded,
                gap_open,
                placed_weights,
                mismatch,
                gap,
            )
            var exact = row_cost(exact_rows, len(rows) - 1, gap_open, placed_weights, mismatch, gap)
            _ = place(rows, len(rows) - 1, centre, forced, band, max_band, gap_open, placed_weights, mismatch, gap)
            var got = row_cost(rows, len(rows) - 1, gap_open, placed_weights, mismatch, gap)
            var floor = 0
            var prices = Costs.affine(mismatch, gap_open, gap)
            for r in range(len(placed) - 1):
                floor += placed_weights[r] * aligner.distance(Span(placed[r]), Span(placed[len(placed) - 1]), prices)
            print(
                "audit step",
                step,
                "sequence",
                order[step],
                "length",
                m,
                "anchored",
                len(chain_j),
                "placed",
                got,
                "exact",
                exact,
                "floor",
                floor,
                file=stderr,
            )
            continue
        var guided = False
        if corridor >= 0 and len(forced) == 0:
            # The new sequence kept inside its co-optimal corridor against its nearest placed relative: at every
            # column boundary, the letters an alignment to that relative within `corridor` of its optimum consumes.
            var guide = -1
            var guide_d = Int.MAX
            for r in range(len(placed) - 1):
                if len(placed[r]) > 0:
                    var d = aligner.distance(Span(placed[r]), Span(placed[len(placed) - 1]))
                    if d < guide_d:
                        guide_d = d
                        guide = r
            if guide >= 0:
                var bounds = _corridor(placed[guide], placed[len(placed) - 1], corridor)
                if len(bounds[0]) > 0:
                    var columns_now = len(rows[0])
                    var j_lo = List[Int](length=columns_now + 1, fill=0)
                    var j_hi = List[Int](length=columns_now + 1, fill=m)
                    var widest = 0
                    var i = 0
                    for c in range(columns_now + 1):
                        j_lo[c] = bounds[0][i]
                        j_hi[c] = bounds[1][i]
                        centre[c] = (j_lo[c] + j_hi[c]) // 2
                        widest = max(widest, (j_hi[c] - j_lo[c] + 1) // 2)
                        if c < columns_now and rows[guide][c] != GAP:
                            i += 1
                    var wide = min(max(band, widest + EDGE_SLACK), max(band, max_band))
                    _ = place(
                        rows,
                        len(rows) - 1,
                        centre,
                        List[Int](),
                        wide,
                        max_band,
                        gap_open,
                        placed_weights,
                        mismatch,
                        gap,
                        j_lo,
                        j_hi,
                    )
                    guided = True
        if guided:
            pass
        elif exact_cells > 0 and len(forced) == 0 and (len(rows[0]) + 1) * (m + 1) <= exact_cells:
            # Exact placement where affordable: the DP unbanded, whatever the anchors say.
            var unbanded = max(len(rows[0]), m) + 1
            _ = place(
                rows,
                len(rows) - 1,
                List[Int](),
                List[Int](),
                unbanded,
                unbanded,
                gap_open,
                placed_weights,
                mismatch,
                gap,
            )
        else:
            _ = place(rows, len(rows) - 1, centre, forced, band, max_band, gap_open, placed_weights, mismatch, gap)
        if interleave > 0 and step >= 2 and step % interleave == 0 and step + 1 < n:
            # A refinement pass over the rows placed so far, before the next sequence is fitted to them.
            rows = refine(
                rows, 1, band, max_band, True, vote_costs, gap_open, False, 0, 0, placed_weights, 1000, mismatch, gap
            )
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


def _two_clusters(seqs: List[List[UInt8]], mut aligner: Aligner) raises -> List[Bool]:
    """The sequences split in two at the top of a UPGMA tree over their edit distances: membership of one side."""
    var n = len(seqs)
    var distance = List[Float64](length=n * n, fill=0.0)
    for a in range(n):
        for b in range(a + 1, n):
            var d = Float64(aligner.distance(Span(seqs[a]), Span(seqs[b])))
            distance[a * n + b] = d
            distance[b * n + a] = d
    var members = List[List[Bool]]()
    var size = List[Int]()
    var alive = List[Bool](length=n, fill=True)
    for a in range(n):
        var m = List[Bool](length=n, fill=False)
        m[a] = True
        members.append(m^)
        size.append(1)
    var left = n
    while left > 2:
        var best = Float64.MAX
        var x = -1
        var y = -1
        for a in range(n):
            if not alive[a]:
                continue
            for b in range(a + 1, n):
                if alive[b] and distance[a * n + b] < best:
                    best = distance[a * n + b]
                    x = a
                    y = b
        for r in range(n):
            if members[y][r]:
                members[x][r] = True
        alive[y] = False
        for r in range(n):
            if alive[r] and r != x:
                var d = (distance[x * n + r] * Float64(size[x]) + distance[y * n + r] * Float64(size[y])) / Float64(
                    size[x] + size[y]
                )
                distance[x * n + r] = d
                distance[r * n + x] = d
        size[x] += size[y]
        left -= 1
    for a in range(n):
        if alive[a]:
            return members[a].copy()
    return List[Bool](length=n, fill=True)


def _merged(
    rows_a: List[List[UInt8]],
    rep_a: Int,
    rows_b: List[List[UInt8]],
    rep_b: Int,
    mut aligner: Aligner,
) raises -> Tuple[List[List[UInt8]], List[List[UInt8]]]:
    """Two MSAs merged into one set of columns, a column of each paired where their representatives' (rows `rep_a`
    and `rep_b`) optimal unit-cost alignment pairs those letters, every other column on its own: the rows of each
    over the joint columns."""
    var col_a = List[Int]()
    var spelled_a = List[UInt8]()
    for c in range(len(rows_a[rep_a])):
        if rows_a[rep_a][c] != GAP:
            col_a.append(c)
            spelled_a.append(rows_a[rep_a][c])
    var col_b = List[Int]()
    var spelled_b = List[UInt8]()
    for c in range(len(rows_b[rep_b])):
        if rows_b[rep_b][c] != GAP:
            col_b.append(c)
            spelled_b.append(rows_b[rep_b][c])
    var found = aligner.align(Span(spelled_a), Span(spelled_b), Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT)
    var runs = cigar_runs(found.cigar)
    # The joint columns: (column of A or -1, column of B or -1), in order.
    var joint_a = List[Int]()
    var joint_b = List[Int]()
    var ca = 0
    var cb = 0
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
                var ta = col_a[i]
                var tb = col_b[j]
                while ca < ta:
                    joint_a.append(ca)
                    joint_b.append(-1)
                    ca += 1
                while cb < tb:
                    joint_a.append(-1)
                    joint_b.append(cb)
                    cb += 1
                joint_a.append(ta)
                joint_b.append(tb)
                ca = ta + 1
                cb = tb + 1
                i += 1
                j += 1
    while ca < len(rows_a[0]):
        joint_a.append(ca)
        joint_b.append(-1)
        ca += 1
    while cb < len(rows_b[0]):
        joint_a.append(-1)
        joint_b.append(cb)
        cb += 1
    var out_a = List[List[UInt8]]()
    for row in rows_a:
        var r = List[UInt8](capacity=len(joint_a))
        for t in range(len(joint_a)):
            r.append(row[joint_a[t]] if joint_a[t] >= 0 else GAP)
        out_a.append(r^)
    var out_b = List[List[UInt8]]()
    for row in rows_b:
        var r = List[UInt8](capacity=len(joint_b))
        for t in range(len(joint_b)):
            r.append(row[joint_b[t]] if joint_b[t] >= 0 else GAP)
        out_b.append(r^)
    return (out_a^, out_b^)


def _tree_rows(seqs: List[List[UInt8]], mut aligner: Aligner) raises -> List[List[UInt8]]:
    """An MSA along the minimum spanning tree of the sequences' edit distances (Prim's, from the longest): each
    sequence joins through its optimal unit-cost alignment to its tree parent, its letters taking the parent's
    columns, the rest columns of their own. Alignments along a tree are always consistent, so every sequence keeps
    its optimal alignment to its nearest relative exactly. Rows in `seqs`' order."""
    var n = len(seqs)
    var longest = 0
    for k in range(1, n):
        if len(seqs[k]) > len(seqs[longest]):
            longest = k
    var near = List[Int](length=n, fill=Int.MAX)
    var parent = List[Int](length=n, fill=-1)
    var done = List[Bool](length=n, fill=False)
    var order = List[Int]()
    var at = longest
    while at >= 0:
        done[at] = True
        order.append(at)
        var next = -1
        for k in range(n):
            if done[k]:
                continue
            var d = aligner.distance(Span(seqs[at]), Span(seqs[k]))
            if d < near[k]:
                near[k] = d
                parent[k] = at
            if next < 0 or near[k] < near[next]:
                next = k
        at = next
    # rows[k] is sequence k's row once it joined; the MSA grows column by column.
    var rows = List[List[UInt8]]()
    for _ in range(n):
        rows.append(List[UInt8]())
    rows[order[0]] = seqs[order[0]].copy()
    var joined = List[Int]()
    joined.append(order[0])
    for t in range(1, len(order)):
        var k = order[t]
        var p = parent[k]
        var placed = List[List[UInt8]]()
        var rep = 0
        for u in range(len(joined)):
            placed.append(rows[joined[u]].copy())
            if joined[u] == p:
                rep = u
        var single = List[List[UInt8]]()
        single.append(seqs[k].copy())
        var both = _merged(placed, rep, single, 0, aligner)
        for u in range(len(joined)):
            rows[joined[u]] = both[0][u].copy()
        rows[k] = both[1][0].copy()
        joined.append(k)
    return rows^


def progressive(
    seqs: List[List[UInt8]],
    passes: Int = 20,
    band: Int = 64,
    max_band: Int = 4096,
    gap_open: Int = 1,
    group_band: Int = 256,
    units: Bool = True,
    even_start: Bool = True,
    one_indel: Int = -2,
    votes: Int = AFFINE_VOTES,
    mismatch: Int = 2,
    gap: Int = 2,
    hard: Int = 0,
    nearest: Bool = False,
    consistency: Bool = False,
    subfamilies: Bool = False,
    shuffle: Int = 0,
    tree_start: Bool = False,
    exact_cells: Int = 0,
    interleave: Int = 0,
    corridor: Int = -1,
) raises -> List[List[UInt8]]:
    """The sequences' MSA (module docstring): one row per sequence, in order, identical sequences sharing their row.
    `passes`, `band`, `max_band`, `gap_open` and `group_band` are the refinement's (see `refine`); `units` builds the
    copy-level start too where the sequences hold a tandem array, and gives the refinement the repeat's period;
    `even_start`, where some sequence repeats, one more start placed and refined with every distinct sequence
    weighing one, then refined weighted (a different path to a local optimum, which is sometimes lower);
    `one_indel`, the seeds' rule inside tandem repeats (`_cascade`); `votes`, the costs of the optimal pairwise
    alignments that vote for placement and centre the refinement's bands: `UNIT_VOTES` (vntr-eval's unit cost),
    `PRICED_VOTES` (gap-affine at the objective's prices: `mismatch`, a gap of k letters `gap_open + gap * k`),
    `AFFINE_VOTES` (vntr-eval's affine cost, 4 and 6 + 2k; the default: on vntr-eval's 149 regions it lowers the
    affine cost over the optima against unit-cost votes, hotspot VNTRs better in 28 and worse in 12, the unit cost
    level) or `TWO_PIECE_VOTES` (minimap2's two-piece costs); `mismatch`
    and `gap`, the prices of a column of two different letters and of a letter against a gap (`refine`); `hard`, how
    many of the best seed sources' anchors hold their letters in placement (`_placed`); `nearest`, sequences placed
    nearest first, not longest first (`_placed`); `consistency`, a letter anchored where most of its votes agree
    (not more than half), each anchor weighing the votes its column won in the chain; `subfamilies`, the subfamily
    start (each of the two clusters at the top of a tree aligned on its own, then merged); `shuffle` > 0, the first
    start's placement in a random order after the longest, seeded by it (a restart); `tree_start`, the start along
    the minimum spanning tree (`_tree_rows`: every sequence optimally aligned to its nearest relative);
    `exact_cells` > 0, each placement exact (unbanded) where its DP has that many cells or fewer; `interleave` > 0,
    a refinement pass over the rows placed so far after every that many placements; `corridor` >= 0, each placement
    kept inside its co-optimal corridor (that many units of slack) against its nearest placed relative."""
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
    var two_piece = Costs.edit()
    if votes == PRICED_VOTES:
        two_piece = Costs.affine(mismatch, gap_open, gap)
    elif votes == AFFINE_VOTES:
        two_piece = Costs.affine(4, 6, 2)
    elif votes == TWO_PIECE_VOTES:
        two_piece = Costs.two_piece(4, 4, 2, 24, 1)
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
    var rows = _placed(
        distinct,
        band,
        max_band,
        gap_open,
        weight,
        one_indel,
        two_piece,
        mismatch,
        gap,
        hard,
        nearest,
        consistency,
        shuffle,
        exact_cells,
        interleave,
        corridor,
    )
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
    # The subfamily start: the two clusters at the top of a tree over the sequences aligned each on its own (a
    # minority whose arrays differ from the majority's keeps its own alignment), merged through their longest
    # members' optimal alignment, then refined together.
    if subfamilies and passes > 0 and len(distinct) >= 6 and len(distinct) <= SUBFAMILY_ROWS:
        var cluster_aligner = Aligner()
        var in_a = _two_clusters(distinct, cluster_aligner)
        var na = 0
        for x in in_a:
            if x:
                na += 1
        if na >= 2 and len(distinct) - na >= 2:
            var seqs_a = List[List[UInt8]]()
            var seqs_b = List[List[UInt8]]()
            var index_a = List[Int]()
            var index_b = List[Int]()
            for d in range(len(distinct)):
                # Each distinct sequence as often as the input holds it, so the clusters keep their weights.
                for _ in range(weight[d]):
                    if in_a[d]:
                        seqs_a.append(distinct[d].copy())
                    else:
                        seqs_b.append(distinct[d].copy())
                if in_a[d]:
                    index_a.append(d)
                else:
                    index_b.append(d)
            var sub_a = progressive(
                seqs_a,
                passes,
                band,
                max_band,
                gap_open,
                group_band,
                units,
                False,
                one_indel,
                votes,
                mismatch,
                gap,
                hard,
                nearest,
                consistency,
                False,
            )
            var sub_b = progressive(
                seqs_b,
                passes,
                band,
                max_band,
                gap_open,
                group_band,
                units,
                False,
                one_indel,
                votes,
                mismatch,
                gap,
                hard,
                nearest,
                consistency,
                False,
            )
            # One row per distinct member: its first occurrence in the cluster's input.
            var rows_a = List[List[UInt8]]()
            var rows_b = List[List[UInt8]]()
            var at = 0
            var rep_a = 0
            for t in range(len(index_a)):
                rows_a.append(sub_a[at].copy())
                if len(distinct[index_a[t]]) > len(distinct[index_a[rep_a]]):
                    rep_a = t
                at += weight[index_a[t]]
            at = 0
            var rep_b = 0
            for t in range(len(index_b)):
                rows_b.append(sub_b[at].copy())
                if len(distinct[index_b[t]]) > len(distinct[index_b[rep_b]]):
                    rep_b = t
                at += weight[index_b[t]]
            var both = _merged(rows_a, rep_a, rows_b, rep_b, cluster_aligner)
            var joint = List[List[UInt8]]()
            for _ in range(len(distinct)):
                joint.append(List[UInt8]())
            for t in range(len(index_a)):
                joint[index_a[t]] = both[0][t].copy()
            for t in range(len(index_b)):
                joint[index_b[t]] = both[1][t].copy()
            joint = refine(
                joint,
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
            if objective(joint, gap_open, weight, mismatch, gap) < objective(rows, gap_open, weight, mismatch, gap):
                rows = joint^
    if tree_start and passes > 0 and len(distinct) > 2 and len(distinct) <= SUBFAMILY_ROWS:
        var tree_aligner = Aligner()
        var tree = _tree_rows(distinct, tree_aligner)
        tree = refine(
            tree,
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
        if objective(tree, gap_open, weight, mismatch, gap) < objective(rows, gap_open, weight, mismatch, gap):
            rows = tree^
    var repeated = len(distinct) < n
    if even_start and repeated and passes > 0 and len(distinct) > 1:
        var even = List[Int](length=len(distinct), fill=1)
        var other = _placed(
            distinct,
            band,
            max_band,
            gap_open,
            even,
            one_indel,
            two_piece,
            mismatch,
            gap,
            hard,
            nearest,
            consistency,
            0,
            exact_cells,
            interleave,
            corridor,
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
                left,
                band,
                max_band,
                gap_open,
                weight,
                one_indel,
                two_piece,
                mismatch,
                gap,
                hard,
                nearest,
                consistency,
                0,
                exact_cells,
                interleave,
                corridor,
            )
            var right_rows = _placed(
                right,
                band,
                max_band,
                gap_open,
                weight,
                one_indel,
                two_piece,
                mismatch,
                gap,
                hard,
                nearest,
                consistency,
                0,
                exact_cells,
                interleave,
                corridor,
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
