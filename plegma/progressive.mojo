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

from karan import Aligner, Costs, Mode, Ties, alignments

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
    # PLEGMA_ORDER (an experiment): `heavy` places the heaviest first (most copies, then longest); `heavy_prim` starts
    # from the heaviest and then takes each time the sequence nearest one already placed.
    var order_text = getenv("PLEGMA_ORDER")
    if order_text == "heavy" or order_text == "heavy_prim":
        var by_weight = List[UInt64]()
        for k in range(n):
            by_weight.append(
                (UInt64((1 << 20) - min(weights[k], (1 << 20) - 1)) << 40)
                | (UInt64((1 << 20) - min(len(seqs[k]), (1 << 20) - 1)) << 20)
                | UInt64(k)
            )
        sort(by_weight)
        order.clear()
        for key in by_weight:
            order.append(Int(key & 0xFFFFF))
    var prim_needed = nearest or order_text == "heavy_prim"
    if prim_needed and n > 2:
        # Prim's order from the first: each next sequence the one nearest any placed one.
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
    # PLEGMA_CENTRE=<k> (an experiment): sequence k placed first and the only voter of every other, a star around it.
    var centre_text = getenv("PLEGMA_CENTRE")
    var centred = centre_text != "" and Int(centre_text) < n
    if centred:
        var centre = Int(centre_text)
        var rest = List[Int]()
        rest.append(centre)
        for k in order:
            if k != centre:
                rest.append(k)
        order = rest^
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
        if centred:
            for r in range(1, len(placed)):
                voters[r] = False
        elif len(placed) > voter_cap:
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


comptime REFERENCE_MEMORY = 2 << 30
"""`on_reference`'s memory cap per pair (karan's default, 80 MB, refuses some divergent pairs of a few kb)."""


def on_reference(
    seqs: List[List[UInt8]], costs: Costs, snap: Bool = False, threads: Int = 1
) raises -> List[List[UInt8]]:
    """The MSA of `seqs` anchored on the first, the reference: every other distinct sequence aligned to it exactly
    (karan, `costs`, gaps shifted left through repeats as variant callers normalize them) and its letters put in the
    reference letters' columns, the letters it inserts in the columns between them (the runs different sequences
    insert at one place aligned to each other the same way, on the longest). So the alignment of the reference with every sequence
    is that pair's own optimal one: what a caller reports against the reference path keeps each event whole, where a
    placement by many partners' votes can split a copy-number change into copy-sized pieces. Rows in `seqs`' order.

    With `threads` above one, the alignments to the reference run over that many threads (karan's batch).

    With `snap`, each insertion or deletion of `SNAP_LENGTH` letters or more then moves, whole, to a place a sequence
    before it (the heavier, as the input repeats them, then the earlier) put an event of its kind and length (within
    a tenth) within `_snap_window` of it, where that costs no more than `_snap_tolerance` mismatches: inside a tandem array the sequences' own optima put
    one copy-number change at places a few copies apart, by their copies' few differences, and the graph then holds
    it as several bubbles whose reads split between them."""
    var n = len(seqs)
    if n == 0:
        return List[List[UInt8]]()
    var distinct = _distinct(seqs)
    var partner = List[List[Int]]()
    var identity = List[Int]()
    for i in range(len(seqs[distinct[0]])):
        identity.append(i)
    partner.append(identity^)
    if threads > 1 and len(distinct) > 2:
        # The alignments are independent: karan's batch spreads them over `threads`, the longest first.
        var references = List[String]()
        var queries = List[String]()
        for d in range(1, len(distinct)):
            references.append(String(unsafe_from_utf8=Span(seqs[distinct[0]])))
            queries.append(String(unsafe_from_utf8=Span(seqs[distinct[d]])))
        var found = alignments(
            references, queries, costs, Mode.GLOBAL, ties=Ties.LEFT, threads=threads, max_memory=REFERENCE_MEMORY
        )
        for d in range(1, len(distinct)):
            partner.append(_pairing(found[d - 1].cigar, len(seqs[distinct[d]])))
    else:
        var aligner = Aligner()
        for d in range(1, len(distinct)):
            partner.append(_onto(aligner, seqs[distinct[0]], seqs[distinct[d]], costs))
    if snap:
        _snap(seqs, distinct, partner)
    var rows_of = _on_reference_columns(seqs, distinct, partner, costs)
    var out = List[List[UInt8]]()
    for k in range(n):
        out.append(rows_of[_index_of(seqs, distinct, k)].copy())
    return out^


comptime SNAP_LENGTH = 10
"""The shortest insertion or deletion `on_reference`'s snapping moves."""


def _snap_window(length: Int) -> Int:
    """How far an event of `length` letters may move."""
    return max(2 * length, 100)


def _snap_tolerance(length: Int) -> Int:
    """The mismatches an event of `length` letters may gain by moving (PLEGMA_SNAP_SLACK=a,b: a + length // b, an
    experiment; 0 by default: only to a place as good for the sequence, so its alignment with the reference stays an
    optimal one and its calls are written as the truth's are)."""
    var slack = getenv("PLEGMA_SNAP_SLACK")
    if slack == "":
        return 0
    var parts = slack.split(",")
    try:
        return Int(String(parts[0])) + length // Int(String(parts[1]))
    except:
        return 0


def _snap(seqs: List[List[UInt8]], distinct: List[Int], mut partner: List[List[Int]]):
    """Each distinct sequence's insertions and deletions (`partner`: each letter's reference letter, -1 inserted)
    moved, whole, to a place an earlier-processed sequence's event of the same kind holds nearby, when the letters
    they pass then pair at no more than the tolerance's extra mismatches; the heavier sequences first."""
    ref c = seqs[distinct[0]]
    var weight = List[Int](length=len(distinct), fill=0)
    for k in range(len(seqs)):
        weight[_index_of(seqs, distinct, k)] += 1
    var order = List[Int]()
    for d in range(1, len(distinct)):
        order.append(d)
    # Heaviest first, the input's order among equals (an insertion sort: a few hundred sequences at most).
    for x in range(1, len(order)):
        var y = x
        while y > 0 and weight[order[y]] > weight[order[y - 1]]:
            var t = order[y]
            order[y] = order[y - 1]
            order[y - 1] = t
            y -= 1
    # The places events were put: an insertion's slot (after reference letter s), a deletion's first letter.
    var insert_slots = List[Int]()
    var insert_lengths = List[Int]()
    var delete_starts = List[Int]()
    var delete_lengths = List[Int]()
    for d in order:
        ref q = seqs[distinct[d]]
        ref p = partner[d]
        # Insertion runs [j0, j1) between reference letters s and s + 1.
        var j = 0
        while j < len(q):
            if p[j] >= 0 or j == 0:
                j += 1
                continue
            var j0 = j
            while j < len(q) and p[j] < 0:
                j += 1
            var j1 = j
            var length = j1 - j0
            var s = p[j0 - 1]
            if s < 0 or j1 >= len(q) or p[j1] != s + 1 or length < SNAP_LENGTH:
                continue
            var best = -1
            var best_delta = _snap_tolerance(length) + 1
            for e in range(len(insert_slots)):
                var slot = insert_slots[e]
                if slot == s or abs(slot - s) > _snap_window(length) or not _same_length(insert_lengths[e], length):
                    continue
                var delta = 0
                var ok = True
                if slot < s:
                    var m = s - slot
                    if j0 - m < 0:
                        continue
                    for i in range(m):
                        if p[j0 - m + i] != slot + 1 + i:
                            ok = False
                            break
                        delta += Int(q[j0 - m + length + i] != c[slot + 1 + i]) - Int(q[j0 - m + i] != c[slot + 1 + i])
                else:
                    var m = slot - s
                    if j1 + m > len(q):
                        continue
                    for i in range(m):
                        if p[j1 + i] != s + 1 + i:
                            ok = False
                            break
                        delta += Int(q[j0 + i] != c[s + 1 + i]) - Int(q[j1 + i] != c[s + 1 + i])
                if ok and delta < best_delta:
                    best_delta = delta
                    best = slot
            if best >= 0:
                if best < s:
                    var m = s - best
                    for i in range(length):
                        p[j0 - m + i] = -1
                    for i in range(m):
                        p[j0 - m + length + i] = best + 1 + i
                else:
                    var m = best - s
                    for i in range(m):
                        p[j0 + i] = s + 1 + i
                    for i in range(length):
                        p[j0 + m + i] = -1
                    # Past the moved run: it is not moved twice.
                    j = j0 + m + length
        # Deletion runs: reference letters [a, a + length) between letters j - 1 and j of the sequence.
        var moved_to = 0
        for k in range(1, len(q)):
            if k <= moved_to or p[k] < 0 or p[k - 1] < 0 or p[k] - p[k - 1] - 1 < SNAP_LENGTH:
                continue
            var a = p[k - 1] + 1
            var length = p[k] - a
            var best = -1
            var best_delta = _snap_tolerance(length) + 1
            for e in range(len(delete_starts)):
                var start = delete_starts[e]
                if start == a or abs(start - a) > _snap_window(length) or not _same_length(delete_lengths[e], length):
                    continue
                var delta = 0
                var ok = True
                if start < a:
                    var m = a - start
                    if k - m < 0:
                        continue
                    for i in range(m):
                        if p[k - m + i] != start + i:
                            ok = False
                            break
                        delta += Int(q[k - m + i] != c[start + length + i]) - Int(q[k - m + i] != c[start + i])
                else:
                    var m = start - a
                    if k + m > len(q) or start + length > len(c):
                        continue
                    for i in range(m):
                        if p[k + i] != a + length + i:
                            ok = False
                            break
                        delta += Int(q[k + i] != c[a + i]) - Int(q[k + i] != c[a + length + i])
                if ok and delta < best_delta:
                    best_delta = delta
                    best = start
            if best >= 0:
                if best < a:
                    var m = a - best
                    for i in range(m):
                        p[k - m + i] = best + length + i
                else:
                    var m = best - a
                    for i in range(m):
                        p[k + i] = a + i
                    moved_to = k + m
        # This sequence's events, as they now stand, are places for the lighter ones.
        var x = 1
        while x < len(q):
            if p[x] < 0 and p[x - 1] >= 0:
                var y = x
                while y < len(q) and p[y] < 0:
                    y += 1
                insert_slots.append(p[x - 1])
                insert_lengths.append(y - x)
                x = y
                continue
            if p[x] >= 0 and p[x - 1] >= 0 and p[x] - p[x - 1] - 1 >= SNAP_LENGTH:
                delete_starts.append(p[x - 1] + 1)
                delete_lengths.append(p[x] - p[x - 1] - 1)
            x += 1


def _same_length(a: Int, b: Int) -> Bool:
    """Two events' lengths within a tenth (5 letters at least): one change, not two changes at one place."""
    return abs(a - b) <= max(5, max(a, b) // 10)


def _distinct(seqs: List[List[UInt8]]) -> List[Int]:
    """The first index of each distinct sequence, in order."""
    var distinct = List[Int]()
    for k in range(len(seqs)):
        var seen = False
        for d in range(len(distinct)):
            if seqs[distinct[d]] == seqs[k]:
                seen = True
                break
        if not seen:
            distinct.append(k)
    return distinct^


def _index_of(seqs: List[List[UInt8]], distinct: List[Int], k: Int) -> Int:
    """Which distinct sequence sequence `k` is."""
    for d in range(len(distinct)):
        if seqs[distinct[d]] == seqs[k]:
            return d
    return -1


def _on_reference_columns(
    seqs: List[List[UInt8]], distinct: List[Int], partner: List[List[Int]], costs: Costs
) raises -> List[List[UInt8]]:
    """The rows of the distinct sequences (`distinct[0]` the reference) whose letters `partner` pairs with reference
    letters (-1: inserted): each letter on its reference letter's column, the runs inserted at one place aligned to
    each other on the longest (`on_reference`, `costs`), so their copies share columns."""
    ref c = seqs[distinct[0]]
    var width = len(c)
    # For an inserted letter, its slot (after reference letter i, -1 for before the first) and its place in that
    # slot's run.
    var slot = List[List[Int]]()
    var place_in = List[List[Int]]()
    var slot_width = List[Int](length=width + 1, fill=0)
    for d in range(len(distinct)):
        var slots = List[Int]()
        var places = List[Int]()
        var last = -1
        var run = 0
        for j in range(len(partner[d])):
            if partner[d][j] >= 0:
                slots.append(-2)
                places.append(0)
                last = partner[d][j]
                run = 0
            else:
                slots.append(last)
                places.append(run)
                run += 1
                slot_width[last + 1] = max(slot_width[last + 1], run)
        slot.append(slots^)
        place_in.append(places^)
    # Where several sequences insert different letters at one place, those runs aligned to each other the same way
    # (on the longest of them), so their copies share columns rather than stand left-justified side by side.
    var start_of = List[Dict[Int, Int]]()
    for d in range(len(distinct)):
        var starts = Dict[Int, Int]()
        for j in range(len(partner[d])):
            if partner[d][j] < 0 and (j == 0 or partner[d][j - 1] >= 0 or slot[d][j - 1] != slot[d][j]):
                starts[slot[d][j] + 1] = j
        start_of.append(starts^)
    for place in range(width + 1):
        if slot_width[place] < 2:
            continue
        var members = List[Int]()
        var pieces = List[List[UInt8]]()
        for d in range(len(distinct)):
            var j0 = start_of[d].get(place)
            if not j0:
                continue
            var j = j0.value()
            var piece = List[UInt8]()
            while j < len(partner[d]) and partner[d][j] < 0 and slot[d][j] + 1 == place:
                piece.append(seqs[distinct[d]][j])
                j += 1
            members.append(d)
            pieces.append(piece^)
        var differ = False
        for k in range(1, len(pieces)):
            if pieces[k] != pieces[0]:
                differ = True
                break
        if not differ:
            continue
        # The longest piece first: the anchor of this place's own alignment.
        var longest = 0
        for k in range(1, len(pieces)):
            if len(pieces[k]) > len(pieces[longest]):
                longest = k
        var ordered = List[List[UInt8]]()
        var at = List[Int]()
        ordered.append(pieces[longest].copy())
        at.append(longest)
        for k in range(len(pieces)):
            if k != longest:
                ordered.append(pieces[k].copy())
                at.append(k)
        var sub = on_reference(ordered, costs)
        slot_width[place] = len(sub[0])
        for t in range(len(at)):
            var d = members[at[t]]
            var j = start_of[d][place]
            for c in range(len(sub[t])):
                if sub[t][c] != GAP:
                    place_in[d][j] = c
                    j += 1
    # Columns: slot before the first reference letter, then each reference letter followed by its slot.
    var column_of_letter = List[Int](length=width, fill=0)
    var slot_start = List[Int](length=width + 1, fill=0)
    var col = 0
    slot_start[0] = 0
    col += slot_width[0]
    for i in range(width):
        column_of_letter[i] = col
        col += 1
        slot_start[i + 1] = col
        col += slot_width[i + 1]
    var rows_of = List[List[UInt8]]()
    for d in range(len(distinct)):
        ref q = seqs[distinct[d]]
        var row = List[UInt8](length=col, fill=GAP)
        for j in range(len(q)):
            if partner[d][j] >= 0:
                row[column_of_letter[partner[d][j]]] = q[j]
            else:
                row[slot_start[slot[d][j] + 1] + place_in[d][j]] = q[j]
        rows_of.append(row^)
    return rows_of^


comptime TREE_K = 15
"""`on_tree`'s k-mers, for each sequence's nearest relative."""


def _kmer_counts(seq: List[UInt8]) -> List[UInt64]:
    """The sequence's k-mers of `TREE_K` ACGT letters as 2-bit codes, sorted, a k-mer as often as it occurs (copies
    of a tandem unit count apart, so two arrays differing in copy number differ here)."""
    var out = List[UInt64]()
    var code = UInt64(0)
    var valid = 0
    var mask = (UInt64(1) << UInt64(2 * TREE_K)) - 1
    for x in seq:
        var b: UInt64
        if x == UInt8(ord("A")):
            b = 0
        elif x == UInt8(ord("C")):
            b = 1
        elif x == UInt8(ord("G")):
            b = 2
        elif x == UInt8(ord("T")):
            b = 3
        else:
            valid = 0
            continue
        code = ((code << 2) | b) & mask
        valid += 1
        if valid >= TREE_K:
            out.append(code)
    sort(out)
    return out^


def _similarity(a: List[UInt64], b: List[UInt64]) -> Float64:
    """The weighted Jaccard index of two sorted k-mer multisets: the k-mers they share (each as often as the rarer
    holds it) over the k-mers either holds."""
    var shared = 0
    var i = 0
    var j = 0
    while i < len(a) and j < len(b):
        if a[i] == b[j]:
            shared += 1
            i += 1
            j += 1
        elif a[i] < b[j]:
            i += 1
        else:
            j += 1
    var union = len(a) + len(b) - shared
    return Float64(shared) / Float64(union) if union > 0 else 0.0


def _onto(mut aligner: Aligner, c: List[UInt8], q: List[UInt8], costs: Costs) raises -> List[Int]:
    """For each letter of `q`, the letter of `c` its exact alignment to `c` pairs it with (-1: inserted)."""
    var found = aligner.align(Span(c), Span(q), costs, Mode.GLOBAL, ties=Ties.LEFT, max_memory=REFERENCE_MEMORY)
    return _pairing(found.cigar, len(q))


def _pairing(cigar: String, q_length: Int) raises -> List[Int]:
    """For each query letter of a global alignment's CIGAR, the reference letter it pairs with (-1: inserted)."""
    var runs = cigar_runs(cigar)
    var out = List[Int](capacity=q_length)
    var i = 0
    for x in range(len(runs.ops)):
        var op = runs.ops[x]
        for _ in range(runs.counts[x]):
            if op == DELETION:
                i += 1
            elif op == INSERTION:
                out.append(-1)
            else:
                out.append(i)
                i += 1
    return out^


def _events(onto: List[Int], parent_on_ref: List[Int], reference_length: Int) -> Int:
    """The insertions and deletions against the reference of a sequence whose letters pair with its parent's letters
    `onto`, the parent's own pairing with the reference `parent_on_ref`: its runs of letters on no reference letter,
    and of reference letters between two consecutive of its own (or past its ends)."""
    var events = 0
    var last = -1
    var inserting = False
    for j in range(len(onto)):
        var r = parent_on_ref[onto[j]] if onto[j] >= 0 else -1
        if r < 0:
            if not inserting:
                events += 1
            inserting = True
            continue
        inserting = False
        if r > last + 1:
            events += 1
        last = r
    if last + 1 < reference_length:
        events += 1
    return events


def on_tree(seqs: List[List[UInt8]], costs: Costs, share: Float64 = 0.0) raises -> List[List[UInt8]]:
    """The MSA of `seqs` along a tree rooted at the first, the reference: each distinct sequence aligned exactly to
    its parent (karan, `costs`, gaps shifted left), its letters taking the reference letters the parent's take, and
    the runs inserted at one place aligned to each other as `on_reference` aligns them. A difference a family of sequences shares is made once, on the
    tree edge into the family, and every member inherits it the same way: a copy-number change most of the panel
    carries is one deletion of the reference for all of them, an allele that adds copies to it that deletion and an
    insertion, where aligning each to the reference on its own (`on_reference`) splits every member's change at
    places of its own.

    The centres, the reference and every sequence at least `share` of the inputs spell (a sequence weighs as many as
    spell it; 0: every one), join first, nearest first (Prim's, on the weighted Jaccard index of k-mer multisets),
    each through its most similar centre already in; then every other through its most similar centre. A family
    hangs one edge below its centre: along a chain of near relatives the small differences of every edge add up to a
    family's change split in pieces. A sequence joins through its centre only where that leaves it fewer insertions
    and deletions against the reference than its own alignment to the reference does (parsimony): where the centre is
    itself a patchwork of copy differences, a member would inherit them and add copies back beside them. Rows in
    `seqs`' order."""
    var n = len(seqs)
    if n == 0:
        return List[List[UInt8]]()
    var distinct = List[Int]()
    var of = List[Int](length=n, fill=-1)
    for k in range(n):
        for d in range(len(distinct)):
            if seqs[distinct[d]] == seqs[k]:
                of[k] = d
                break
        if of[k] < 0:
            of[k] = len(distinct)
            distinct.append(k)
    var m = len(distinct)
    var weight = List[Int](length=m, fill=0)
    for k in range(n):
        weight[of[k]] += 1
    var centre = List[Bool](length=m, fill=False)
    for d in range(m):
        centre[d] = d == 0 or Float64(weight[d]) >= share * Float64(n)
    var kmers = List[List[UInt64]]()
    for d in range(m):
        kmers.append(_kmer_counts(seqs[distinct[d]]))
    # The centres heaviest first from the reference (the commonest allele takes its own alignment to the reference,
    # a lighter one its relation to a heavier), each joining through the most similar centre already in; the others
    # after, each through its most similar centre. With every sequence a centre (`share` 0), Prim's order instead.
    var parent = List[Int](length=m, fill=0)
    var near = List[Float64](length=m, fill=-1.0)
    var done = List[Bool](length=m, fill=False)
    var order = List[Int]()
    var at = 0
    while at >= 0:
        done[at] = True
        order.append(at)
        var next = -1
        for d in range(m):
            if done[d]:
                continue
            var s = _similarity(kmers[at], kmers[d])
            if s > near[d]:
                near[d] = s
                parent[d] = at
            if centre[d] and (
                next < 0
                or (share > 0.0 and weight[d] > weight[next])
                or ((share == 0.0 or weight[d] == weight[next]) and near[d] > near[next])
            ):
                next = d
        at = next
    for d in range(m):
        if not done[d]:
            order.append(d)
    # Each sequence's letters' reference letters (-1: inserted against the reference), through its parents.
    var on_ref = List[List[Int]]()
    for _ in range(m):
        on_ref.append(List[Int]())
    for i in range(len(seqs[distinct[0]])):
        on_ref[0].append(i)
    var reference_length = len(seqs[distinct[0]])
    var aligner = Aligner()
    for t in range(1, m):
        var d = order[t]
        ref q = seqs[distinct[d]]
        # Through its most similar centre or straight onto the reference, whichever makes it differ from the
        # reference by fewer events (the reference on a tie): a family's shared change is inherited where that is
        # the simpler account, and an allele whose centre's own changes it would inherit and undo stays apart.
        var p = parent[d]
        var onto = _onto(aligner, seqs[distinct[p]], q, costs)
        if p != 0:
            var direct = _onto(aligner, seqs[distinct[0]], q, costs)
            if _events(direct, on_ref[0], reference_length) <= _events(onto, on_ref[p], reference_length):
                p = 0
                onto = direct^
        for j in range(len(q)):
            on_ref[d].append(on_ref[p][onto[j]] if onto[j] >= 0 else -1)
    var rows_of = _on_reference_columns(seqs, distinct, on_ref, costs)
    var out = List[List[UInt8]]()
    for k in range(n):
        out.append(rows_of[of[k]].copy())
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
