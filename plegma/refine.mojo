"""An MSA refined on its sum of pairs: each distinct sequence realigned in turn to the profile of the others.

The sum of pairs at unit costs, over every pair of distinct sequences, counts the columns where the two rows
differ, a column where both are gaps costing nothing: what vntr-eval's all_cost_over_opt divides by the pairs'
optimal costs. With the other rows fixed, a sequence's best row is exact by a DP over the columns: its letter
against column c costs the others that do not hold that letter there (their gaps included), its gap there costs
the others that hold a letter, and its letter in a column of its own costs every other row. The DP keeps within
`band` of the sequence's present place, which is one of its choices, so no pass raises the sum. Passes go on until
one gains nothing. On a hotspot VNTR whose POA (abPOA's, or plegma's longest first) cost 1.29 times the optima,
three passes took it to 1.13.
"""

from karan import Aligner, Costs, Mode, Ties

from .align import anchored_alignment
from .backbone import Level, MUMS
from .chain import ChainCosts
from .graph import DELETION, INSERTION, cigar_runs


comptime GAP = UInt8(ord("-"))
comptime SYMBOLS = 6
"""A, C, G, T, any other letter, and the gap."""
comptime BIG = Int32(1 << 29)
comptime MATCHED = UInt8(0)
comptime GAPPED = UInt8(1)
comptime OWN = UInt8(2)
comptime EDGE = 4


@always_inline
def _symbol(letter: UInt8) -> Int:
    var lower = letter | 0x20
    if letter == GAP:
        return 5
    if lower == UInt8(ord("a")):
        return 0
    if lower == UInt8(ord("c")):
        return 1
    if lower == UInt8(ord("g")):
        return 2
    if lower == UInt8(ord("t")):
        return 3
    return 4


@always_inline
def _w(weights: List[Int], r: Int) -> Int:
    """Row r's weight: how many sequences it stands for (1 where `weights` is empty)."""
    return weights[r] if len(weights) > 0 else 1


def gap_opens(rows: List[List[UInt8]], weights: List[Int] = List[Int]()) -> Int:
    """The gaps the rows' pairs open, summed over every pair: in each pair's alignment (columns where both are gaps
    left out), the runs of columns where just one of the two has a gap; a pair counts the product of the two rows'
    `weights` (1 each where empty)."""
    var total = 0
    var n = len(rows)
    for a in range(n):
        for b in range(a + 1, n):
            var w = _w(weights, a) * _w(weights, b)
            var state = 0
            for c in range(len(rows[a])):
                var x = rows[a][c] == GAP
                var y = rows[b][c] == GAP
                if x and y:
                    continue
                var now = 1 if x else (2 if y else 0)
                if now != 0 and now != state:
                    total += w
                state = now
    return total


def row_cost(rows: List[List[UInt8]], s: Int, gap_open: Int, weights: List[Int] = List[Int]()) -> Int:
    """Row `s`'s exact share of the objective: over every other row, the columns where the two differ (both gaps
    cost nothing), plus `gap_open` for every gap the pair opens, times the other row's weight (`weights`, 1 each
    where empty)."""
    var total = 0
    for r in range(len(rows)):
        if r == s:
            continue
        var pair = 0
        var state = 0
        for c in range(len(rows[s])):
            var x = rows[s][c]
            var y = rows[r][c]
            var gx = x == GAP
            var gy = y == GAP
            if gx and gy:
                continue
            if x != y:
                pair += 1
            var now = 1 if gx else (2 if gy else 0)
            if now != 0 and now != state:
                pair += gap_open
            state = now
        total += pair * _w(weights, r)
    return total


def sum_of_pairs(rows: List[List[UInt8]], weights: List[Int] = List[Int]()) -> Int:
    """The rows' sum of pairs at unit costs (see the module docstring), a pair counting the product of the two rows'
    `weights` (1 each where empty). Letters other than ACGT count as one."""
    if len(rows) == 0:
        return 0
    var n = len(rows)
    var weight = 0
    var squares = 0
    for r in range(n):
        weight += _w(weights, r)
        squares += _w(weights, r) * _w(weights, r)
    var pairs = (weight * weight - squares) // 2
    var total = 0
    var counts = List[Int](length=SYMBOLS, fill=0)
    var counts_sq = List[Int](length=SYMBOLS, fill=0)
    for c in range(len(rows[0])):
        for x in range(SYMBOLS):
            counts[x] = 0
            counts_sq[x] = 0
        for r in range(n):
            var w = _w(weights, r)
            counts[_symbol(rows[r][c])] += w
            counts_sq[_symbol(rows[r][c])] += w * w
        var same = 0
        for x in range(SYMBOLS):
            same += (counts[x] * counts[x] - counts_sq[x]) // 2
        total += pairs - same
    return total


comptime PARTNERS = 24
"""Past this many other distinct sequences, each sequence is aligned to its nearest this many alone (by unit edit
distance), as are those it is nearest to: the anchors' votes need a few dozen, and all pairs cost the square."""


def _pair_runs(
    mut aligner: Aligner, partner: ImmSpan[UInt8, _], query: ImmSpan[UInt8, _], costs: Costs, seeded: Bool
) raises -> Tuple[List[UInt8], List[Int]]:
    """The runs of `query`'s alignment to `partner`: optimal under `costs`, or with `seeded` through anchors."""
    if seeded:
        var levels = List[Level]()
        for k in [51, 31, 21]:
            levels.append(Level(MUMS, k))
        var through = anchored_alignment(
            aligner,
            partner,
            query,
            levels,
            ChainCosts(substitution=0.1),
            costs,
            0,
            True,
            repeat_k=12,
            repeat_distance=500,
            mask_depth=1,
        )
        return (through.ops.copy(), through.counts.copy())
    var found = aligner.align(partner, query, costs, Mode.GLOBAL, ties=Ties.LEFT)
    var runs = cigar_runs(found.cigar)
    return (runs.ops.copy(), runs.counts.copy())


struct Partners(Movable):
    """Each distinct sequence's alignment to its partners (every other, or past PARTNERS others its nearest), as
    the letter of the partner each of its letters faces: `facing[s][r][j]` (-1 where letter j is inserted; the list
    empty where r is no partner of s). They do not depend on the MSA, so they are made once."""

    var facing: List[List[List[Int32]]]

    def __init__(out self, spelled: List[List[UInt8]], costs: Costs = Costs.edit(), seeded: Bool = False) raises:
        """Optimal alignments under `costs` (karan, indels placed left); with `seeded`, through anchors instead
        (`anchored_alignment`): MUMs of 51, 31 and 21 letters outside tandem repeats pin the flanks, then inside each
        gap between them, a repeat array, MUMs unique within the gap's two pieces, which sit on copies that carry
        something of their own and so pair homologous copies; only the stretches between anchors are aligned for
        least cost."""
        var n = len(spelled)
        self.facing = List[List[List[Int32]]]()
        for _ in range(n):
            var row = List[List[Int32]]()
            for _ in range(n):
                row.append(List[Int32]())
            self.facing.append(row^)
        var aligner = Aligner()
        # Which pairs to align: all, or past PARTNERS each sequence's nearest.
        var wanted = List[Bool](length=n * n, fill=n - 1 <= PARTNERS)
        if n - 1 > PARTNERS:
            var distance = List[Int](length=n * n, fill=0)
            for a in range(n):
                for b in range(a + 1, n):
                    var d = aligner.distance(Span(spelled[a]), Span(spelled[b]))
                    distance[a * n + b] = d
                    distance[b * n + a] = d
            var keys = List[UInt64]()
            for a in range(n):
                keys.clear()
                for b in range(n):
                    if b != a:
                        keys.append((UInt64(distance[a * n + b]) << 32) | UInt64(b))
                sort(keys)
                for t in range(PARTNERS):
                    var b = Int(keys[t] & 0xFFFFFFFF)
                    wanted[a * n + b] = True
                    wanted[b * n + a] = True
        for a in range(n):
            for b in range(a + 1, n):
                if not wanted[a * n + b] or len(spelled[a]) == 0 or len(spelled[b]) == 0:
                    continue
                var ab = List[Int32](length=len(spelled[a]), fill=-1)
                var ba = List[Int32](length=len(spelled[b]), fill=-1)
                var runs = _pair_runs(aligner, Span(spelled[b]), Span(spelled[a]), costs, seeded)
                var i = 0
                var j = 0
                for r in range(len(runs[0])):
                    var op = runs[0][r]
                    for _ in range(runs[1][r]):
                        if op == DELETION:
                            i += 1
                        elif op == INSERTION:
                            j += 1
                        else:
                            ab[j] = Int32(i)
                            ba[i] = Int32(j)
                            i += 1
                            j += 1
                self.facing[a][b] = ab^
                self.facing[b][a] = ba^


def _anchored_centre(rows: List[List[UInt8]], s: Int, partners: Partners, mut centre: List[Int]):
    """Where each column's band is centred for row `s`: the sequence's letters before the column as its confident
    anchors put them. Each other row votes, through its optimal pairwise alignment to the sequence, a column for
    every letter of the sequence; a letter is anchored where more than half of the rows vote one column, and the
    anchors kept are the most that rise in both letter and column (a longest increasing chain). Between them the
    centre is interpolated."""
    var columns = len(rows[0])
    var n_rows = len(rows)
    # Each row's letters' columns.
    var column_of = List[List[Int32]]()
    var n = 0
    for r in range(n_rows):
        var cols = List[Int32]()
        for c in range(columns):
            if rows[r][c] != GAP:
                cols.append(Int32(c))
        if r == s:
            n = len(cols)
        column_of.append(cols^)
    var anchor_j = List[Int]()
    var anchor_c = List[Int]()
    var votes = List[Int32]()
    var voters = 0
    for r in range(n_rows):
        if r != s and len(partners.facing[s][r]) > 0:
            voters += 1
    for j in range(n):
        votes.clear()
        for r in range(n_rows):
            if r == s or len(partners.facing[s][r]) == 0:
                continue
            var k = Int(partners.facing[s][r][j])
            if k >= 0:
                votes.append(column_of[r][k])
        if len(votes) * 2 <= voters:
            continue
        sort(votes)
        var best = -1
        var best_count = 0
        var run = 0
        for t in range(len(votes)):
            run = run + 1 if t > 0 and votes[t] == votes[t - 1] else 1
            if run > best_count:
                best_count = run
                best = Int(votes[t])
        if best_count * 2 > voters:
            anchor_j.append(j)
            anchor_c.append(best)
    # The longest chain rising in column (letters already rise): patience sorting.
    var tails = List[Int]()
    var back = List[Int](length=len(anchor_c), fill=-1)
    for t in range(len(anchor_c)):
        var lo = 0
        var hi = len(tails)
        while lo < hi:
            var mid = (lo + hi) // 2
            if anchor_c[tails[mid]] < anchor_c[t]:
                lo = mid + 1
            else:
                hi = mid
        if lo > 0:
            back[t] = tails[lo - 1]
        if lo == len(tails):
            tails.append(t)
        else:
            tails[lo] = t
    var chain = List[Int]()
    var t = tails[len(tails) - 1] if len(tails) > 0 else -1
    while t >= 0:
        chain.append(t)
        t = back[t]
    # Points (column, letters before it): the start, each anchor's column and the one after it, the end.
    var point_c = List[Int]()
    var point_j = List[Int]()
    point_c.append(0)
    point_j.append(0)
    for q in reversed(range(len(chain))):
        var a = chain[q]
        if anchor_c[a] > point_c[len(point_c) - 1] and anchor_j[a] >= point_j[len(point_j) - 1]:
            point_c.append(anchor_c[a])
            point_j.append(anchor_j[a])
            point_c.append(anchor_c[a] + 1)
            point_j.append(anchor_j[a] + 1)
    if point_c[len(point_c) - 1] < columns:
        point_c.append(columns)
        point_j.append(n)
    centre.resize(columns + 1, 0)
    var p = 0
    for c in range(columns + 1):
        while p + 1 < len(point_c) - 1 and point_c[p + 1] <= c:
            p += 1
        var c0 = point_c[p]
        var c1 = point_c[min(p + 1, len(point_c) - 1)]
        var j0 = point_j[p]
        var j1 = point_j[min(p + 1, len(point_j) - 1)]
        if c1 <= c0:
            centre[c] = min(n, max(0, j0))
        else:
            centre[c] = min(n, max(0, j0 + (c - c0) * (j1 - j0) // (c1 - c0)))


def _realign(
    mut rows: List[List[UInt8]],
    s: Int,
    band: Int,
    max_band: Int,
    centre: List[Int],
    gap_open: Int32,
    mut cost: List[Int32],
    mut moves: List[UInt8],
    forced: List[Int] = List[Int](),
    always: Bool = False,
    weights: List[Int] = List[Int](),
) -> Tuple[Int, Int]:
    """Row `s` realigned to the others' profile within `band` letters of its place (or of `centre`, when given),
    widened while the best path touches the band's edge, up to `max_band`; the rows rewritten (all-gap columns
    dropped), unless the new row would cost more than the present one.

    Three states per cell: the sequence's letter in a column (M), its gap there (G), its letter in a column of its
    own (O). Besides the unit costs, every gap a pair newly opens costs `gap_open`, counted from the columns' profile:
    the sequence's gap opening against the rows holding a letter there, a row's own gap starting against the
    sequence's letter, the sequence's own column against every row holding a letter beside it. Those counts only
    approximate the pairs' own gaps (a pair skips the columns where both are gaps), so the DP proposes, and the
    new row is kept only if its exact cost (`row_cost`) is no more than the present one's (or always, with `always`).
    `forced[j]`, where given and not -1, is the column letter j must take. Every other row counts its weight
    (`weights`, 1 each where empty). Returns the exact costs before and after (equal when the row stays)."""
    var has_forced = len(forced) > 0
    var columns = len(rows[0])
    var others = Int32(0)
    for r in range(len(rows)):
        if r != s:
            others += Int32(_w(weights, r))
    # The others' symbol counts per column, and per column how many of them hold a letter, start a gap (a letter
    # in the column before), or end one (a gap before).
    var counts = List[Int32](length=columns * SYMBOLS, fill=0)
    var held = List[Int32](length=columns, fill=0)
    var starts = List[Int32](length=columns, fill=0)
    var ends = List[Int32](length=columns, fill=0)
    for r in range(len(rows)):
        if r == s:
            continue
        var was_letter = True
        var w = Int32(_w(weights, r))
        for c in range(columns):
            var x = rows[r][c]
            counts[c * SYMBOLS + _symbol(x)] += w
            if x != GAP:
                held[c] += w
                if not was_letter:
                    ends[c] += w
            elif was_letter:
                starts[c] += w
            was_letter = x != GAP
    var letters = List[UInt8]()
    var before = List[Int](length=columns + 1, fill=0)
    for c in range(columns):
        before[c] = len(letters)
        if rows[s][c] != GAP:
            letters.append(rows[s][c])
    before[columns] = len(letters)
    var n = len(letters)
    var mid = centre.copy() if len(centre) == columns + 1 else before.copy()
    var g = gap_open
    var b = band
    var steps = List[Int]()
    var mine = List[UInt8]()
    while True:
        var width = 2 * b + 1
        var cells = (columns + 1) * width
        cost.resize(3 * cells, BIG)
        moves.resize(3 * cells, MATCHED)
        # State t's cell (c, j) sits at t * cells + c * width + j - (mid[c] - b).
        for c in range(columns + 1):
            var lo = mid[c] - b
            var base = c * width
            var x = c - 1
            var lo_prev = mid[c - 1] - b if c > 0 else 0
            var prev = (c - 1) * width
            # Opening against the rows holding a letter in this column (before the first: every row).
            var held_here = held[x] if c > 0 else others
            for k in range(width):
                var j = lo + k
                var m_best = BIG
                var m_how = MATCHED
                var g_best = BIG
                var g_how = MATCHED
                var o_best = BIG
                var o_how = MATCHED
                if j >= 0 and j <= n:
                    if c == 0 and j == 0:
                        m_best = 0
                    if c > 0:
                        # M: letter j - 1 in column x.
                        if j > 0 and (not has_forced or forced[j - 1] < 0 or forced[j - 1] == x):
                            var km = j - 1 - lo_prev
                            if km >= 0 and km < width:
                                var unit = others - counts[x * SYMBOLS + _symbol(letters[j - 1])]
                                var from_m = cost[prev + km]
                                var from_g = cost[cells + prev + km]
                                var from_o = cost[2 * cells + prev + km]
                                if from_m < BIG and from_m + g * starts[x] < m_best:
                                    m_best = from_m + g * starts[x]
                                    m_how = MATCHED
                                if from_g < BIG and from_g + g * starts[x] < m_best:
                                    m_best = from_g + g * starts[x]
                                    m_how = GAPPED
                                if from_o < BIG and from_o < m_best:
                                    m_best = from_o
                                    m_how = OWN
                                if m_best < BIG:
                                    m_best += unit
                        # G: the sequence gapped in column x.
                        var kp = j - lo_prev
                        if kp >= 0 and kp < width:
                            var from_m = cost[prev + kp]
                            var from_g = cost[cells + prev + kp]
                            var from_o = cost[2 * cells + prev + kp]
                            if from_m < BIG and from_m + g * held[x] < g_best:
                                g_best = from_m + g * held[x]
                                g_how = MATCHED
                            if from_g < BIG and from_g + g * ends[x] < g_best:
                                g_best = from_g + g * ends[x]
                                g_how = GAPPED
                            if from_o < BIG and from_o + g * held[x] < g_best:
                                g_best = from_o + g * held[x]
                                g_how = OWN
                            if g_best < BIG:
                                g_best += held[x]
                    # O: letter j - 1 in a column of its own after column x.
                    if j > 0 and k > 0 and (not has_forced or forced[j - 1] < 0):
                        var from_m = cost[base + k - 1]
                        var from_g = cost[cells + base + k - 1]
                        var from_o = cost[2 * cells + base + k - 1]
                        if from_m < BIG and from_m + g * held_here < o_best:
                            o_best = from_m + g * held_here
                            o_how = MATCHED
                        if from_g < BIG and from_g + g * held_here < o_best:
                            o_best = from_g + g * held_here
                            o_how = GAPPED
                        if from_o < BIG and from_o < o_best:
                            o_best = from_o
                            o_how = OWN
                        if o_best < BIG:
                            o_best += others
                cost[base + k] = m_best
                moves[base + k] = m_how
                cost[cells + base + k] = g_best
                moves[cells + base + k] = g_how
                cost[2 * cells + base + k] = o_best
                moves[2 * cells + base + k] = o_how
        var end = n - (mid[columns] - b)
        var final_state = 0
        var final_cost = BIG
        if end >= 0 and end < width:
            for t in range(3):
                if cost[t * cells + columns * width + end] < final_cost:
                    final_cost = cost[t * cells + columns * width + end]
                    final_state = t
        if final_cost >= BIG:
            if b >= max_band or b >= max(n, columns):
                var stays = row_cost(rows, s, Int(gap_open), weights)
                return (stays, stays)
            b = min(2 * b, max_band)
            continue
        steps.clear()
        mine.clear()
        var tight = False
        var c = columns
        var j = n
        var state = final_state
        while c > 0 or j > 0:
            var k = j - (mid[c] - b)
            if (k < EDGE and mid[c] - b > 0) or (width - 1 - k < EDGE and mid[c] + b < n):
                tight = True
            var came = Int(moves[state * cells + c * width + k])
            if state == 0:
                steps.append(c - 1)
                mine.append(letters[j - 1])
                c -= 1
                j -= 1
            elif state == 1:
                steps.append(c - 1)
                mine.append(GAP)
                c -= 1
            else:
                steps.append(-1)
                mine.append(letters[j - 1])
                j -= 1
            state = came
        if tight and b < max_band and b < max(n, columns):
            b = min(2 * b, max_band)
            continue
        break
    var fresh = List[List[UInt8]]()
    for r in range(len(rows)):
        var row = List[UInt8](capacity=len(steps))
        for t in reversed(range(len(steps))):
            var old = steps[t]
            if r == s:
                row.append(mine[t])
            elif old >= 0:
                row.append(rows[r][old])
            else:
                row.append(GAP)
        fresh.append(row^)
    var kept = List[Int]()
    for t in range(len(steps)):
        for r in range(len(fresh)):
            if fresh[r][t] != GAP:
                kept.append(t)
                break
    if len(kept) < len(steps):
        for r in range(len(fresh)):
            var row = List[UInt8](capacity=len(kept))
            for t in kept:
                row.append(fresh[r][t])
            fresh[r] = row^
    var before_cost = row_cost(rows, s, Int(gap_open), weights)
    var after_cost = row_cost(fresh, s, Int(gap_open), weights)
    if after_cost > before_cost and not always:
        return (before_cost, before_cost)
    rows = fresh^
    return (before_cost, after_cost)


def _objective(rows: List[List[UInt8]], gap_open: Int, weights: List[Int] = List[Int]()) -> Int:
    if gap_open > 0:
        return sum_of_pairs(rows, weights) + gap_open * gap_opens(rows, weights)
    return sum_of_pairs(rows, weights)


def _tree_groups(rows: List[List[UInt8]]) -> List[List[Bool]]:
    """The clusters of a UPGMA tree over the rows (distance: the columns where two rows differ, both gaps left out),
    each as membership: every merge but the last, so every split of the tree with two or more rows on one side."""
    var n = len(rows)
    var groups = List[List[Bool]]()
    if n < 3:
        return groups^
    var distance = List[Float64](length=n * n, fill=0.0)
    for a in range(n):
        for b in range(a + 1, n):
            var d = 0
            for c in range(len(rows[a])):
                if rows[a][c] != rows[b][c]:
                    d += 1
            distance[a * n + b] = Float64(d)
            distance[b * n + a] = Float64(d)
    var members = List[List[Bool]]()
    var size = List[Int]()
    var alive = List[Bool](length=n, fill=True)
    for a in range(n):
        var m = List[Bool](length=n, fill=False)
        m[a] = True
        members.append(m^)
        size.append(1)
    for _ in range(n - 2):
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
        groups.append(members[x].copy())
    return groups^


def _length_groups(rows: List[List[UInt8]], period: Int) -> List[List[Bool]]:
    """Copy-number splits: the rows whose sequences are longer than each step of half a `period` or more between
    their sorted lengths, against the rest. Sequences of a region that differ by whole copies differ in length by
    multiples of the period, and a gap of whole copies is shared by every shorter sequence alike."""
    var n = len(rows)
    var length = List[Int](length=n, fill=0)
    for r in range(n):
        for x in rows[r]:
            if x != GAP:
                length[r] += 1
    var sorted_lengths = length.copy()
    sort(sorted_lengths)
    var groups = List[List[Bool]]()
    for t in range(1, n):
        if (sorted_lengths[t] - sorted_lengths[t - 1]) * 2 < period:
            continue
        var group = List[Bool](length=n, fill=False)
        for r in range(n):
            group[r] = length[r] >= sorted_lengths[t]
        groups.append(group^)
    return groups^


def _merge_groups(
    mut rows: List[List[UInt8]],
    in_a: List[Bool],
    band: Int,
    max_band: Int,
    gap_open: Int32,
    mut cost: List[Int32],
    mut moves: List[UInt8],
    weights: List[Int] = List[Int](),
) -> Tuple[Int, Int]:
    """The rows in `in_a` (A) and the others (B) realigned as two profiles, each keeping its own columns' order and
    contents, within `band` columns of their present merge (widened up to `max_band` while the path touches the
    band's edge). Pairing A's column with B's costs the pairs across that differ there; A's column alone costs its
    letters times B's rows, and B's alike: at unit costs the DP is exact, since the pairs within A or within B do not
    change. A gap run opening across costs `gap_open` for each pair of a letter there and a letter in the other
    side's column before, an approximation; the new rows are kept only if the exact objective (`_objective`) falls.
    Every row counts its weight (`weights`, 1 each where empty). Returns the objective before and after (equal when
    the rows stay)."""
    var columns = len(rows[0])
    var n_a = 0
    var n_b = 0
    for r in range(len(rows)):
        if in_a[r]:
            n_a += _w(weights, r)
        else:
            n_b += _w(weights, r)
    var before_cost = _objective(rows, Int(gap_open), weights)
    if n_a == 0 or n_b == 0:
        return (before_cost, before_cost)
    # Each side's columns (those where one of its rows holds a letter), their symbol counts and letters held, and
    # where each A column sits in the present merge: the B columns before it.
    var a_cols = List[Int]()
    var b_cols = List[Int]()
    var mid = List[Int]()
    for c in range(columns):
        var has_a = False
        var has_b = False
        for r in range(len(rows)):
            if rows[r][c] != GAP:
                if in_a[r]:
                    has_a = True
                else:
                    has_b = True
        if has_a:
            mid.append(len(b_cols))
            a_cols.append(c)
        if has_b:
            b_cols.append(c)
    var la = len(a_cols)
    var lb = len(b_cols)
    mid.append(lb)
    var counts_a = List[Int32](length=la * SYMBOLS, fill=0)
    var counts_b = List[Int32](length=lb * SYMBOLS, fill=0)
    var held_a = List[Int32](length=la, fill=0)
    var held_b = List[Int32](length=lb, fill=0)
    for r in range(len(rows)):
        var w = Int32(_w(weights, r))
        if in_a[r]:
            for i in range(la):
                var x = rows[r][a_cols[i]]
                counts_a[i * SYMBOLS + _symbol(x)] += w
                if x != GAP:
                    held_a[i] += w
        else:
            for j in range(lb):
                var x = rows[r][b_cols[j]]
                counts_b[j * SYMBOLS + _symbol(x)] += w
                if x != GAP:
                    held_b[j] += w
    var na = Int32(n_a)
    var nb = Int32(n_b)
    var g = gap_open
    # The costs must fit: every pair across, every column, and its openings.
    if n_a * n_b * (1 + 2 * Int(g)) * (la + lb + 1) >= Int(BIG):
        return (before_cost, before_cost)
    var b = band
    var path_a = List[Int]()
    var path_b = List[Int]()
    while True:
        var width = 2 * b + 1
        var cells = (la + 1) * width
        cost.resize(3 * cells, BIG)
        moves.resize(3 * cells, MATCHED)
        # State t's cell (i, j) sits at t * cells + i * width + j - (mid[i] - b): M pairs A's column i - 1 with B's
        # j - 1, A (GAPPED) takes A's column i - 1 alone, B (OWN) takes B's column j - 1 alone.
        for i in range(la + 1):
            var lo = mid[i] - b
            var base = i * width
            var lo_prev = mid[i - 1] - b if i > 0 else 0
            var prev = (i - 1) * width
            var a_before = held_a[i - 1] if i > 0 else na
            for k in range(width):
                var j = lo + k
                var m_best = BIG
                var m_how = MATCHED
                var g_best = BIG
                var g_how = MATCHED
                var o_best = BIG
                var o_how = MATCHED
                if j >= 0 and j <= lb:
                    if i == 0 and j == 0:
                        m_best = 0
                    var b_before = held_b[j - 1] if j > 0 else nb
                    if i > 0:
                        if j > 0:
                            var km = j - 1 - lo_prev
                            if km >= 0 and km < width:
                                var same = Int32(0)
                                for x in range(SYMBOLS):
                                    same += counts_a[(i - 1) * SYMBOLS + x] * counts_b[(j - 1) * SYMBOLS + x]
                                var pair = na * nb - same
                                var from_m = cost[prev + km]
                                var from_g = cost[cells + prev + km]
                                var from_o = cost[2 * cells + prev + km]
                                if from_m < m_best:
                                    m_best = from_m
                                    m_how = MATCHED
                                if from_g < m_best:
                                    m_best = from_g
                                    m_how = GAPPED
                                if from_o < m_best:
                                    m_best = from_o
                                    m_how = OWN
                                if m_best < BIG:
                                    m_best += pair
                        var kp = j - lo_prev
                        if kp >= 0 and kp < width:
                            var opening = g * held_a[i - 1] * b_before
                            var from_m = cost[prev + kp]
                            var from_g = cost[cells + prev + kp]
                            var from_o = cost[2 * cells + prev + kp]
                            if from_m < BIG and from_m + opening < g_best:
                                g_best = from_m + opening
                                g_how = MATCHED
                            if from_g < g_best:
                                g_best = from_g
                                g_how = GAPPED
                            if from_o < BIG and from_o + opening < g_best:
                                g_best = from_o + opening
                                g_how = OWN
                            if g_best < BIG:
                                g_best += held_a[i - 1] * nb
                    if j > 0 and k > 0:
                        var opening = g * held_b[j - 1] * a_before
                        var from_m = cost[base + k - 1]
                        var from_g = cost[cells + base + k - 1]
                        var from_o = cost[2 * cells + base + k - 1]
                        if from_m < BIG and from_m + opening < o_best:
                            o_best = from_m + opening
                            o_how = MATCHED
                        if from_g < BIG and from_g + opening < o_best:
                            o_best = from_g + opening
                            o_how = GAPPED
                        if from_o < o_best:
                            o_best = from_o
                            o_how = OWN
                        if o_best < BIG:
                            o_best += held_b[j - 1] * na
                cost[base + k] = m_best
                moves[base + k] = m_how
                cost[cells + base + k] = g_best
                moves[cells + base + k] = g_how
                cost[2 * cells + base + k] = o_best
                moves[2 * cells + base + k] = o_how
        var end = lb - (mid[la] - b)
        var final_state = 0
        var final_cost = BIG
        if end >= 0 and end < width:
            for t in range(3):
                if cost[t * cells + la * width + end] < final_cost:
                    final_cost = cost[t * cells + la * width + end]
                    final_state = t
        if final_cost >= BIG:
            if b >= max_band or b >= max(la, lb):
                return (before_cost, before_cost)
            b = min(2 * b, max_band)
            continue
        path_a.clear()
        path_b.clear()
        var tight = False
        var i = la
        var j = lb
        var state = final_state
        while i > 0 or j > 0:
            var k = j - (mid[i] - b)
            if (k < EDGE and mid[i] - b > 0) or (width - 1 - k < EDGE and mid[i] + b < lb):
                tight = True
            var came = Int(moves[state * cells + i * width + k])
            if state == 0:
                path_a.append(a_cols[i - 1])
                path_b.append(b_cols[j - 1])
                i -= 1
                j -= 1
            elif state == 1:
                path_a.append(a_cols[i - 1])
                path_b.append(-1)
                i -= 1
            else:
                path_a.append(-1)
                path_b.append(b_cols[j - 1])
                j -= 1
            state = came
        if tight and b < max_band and b < max(la, lb):
            b = min(2 * b, max_band)
            continue
        break
    var fresh = List[List[UInt8]]()
    for r in range(len(rows)):
        var row = List[UInt8](capacity=len(path_a))
        for t in reversed(range(len(path_a))):
            var c = path_a[t] if in_a[r] else path_b[t]
            row.append(rows[r][c] if c >= 0 else GAP)
        fresh.append(row^)
    var after_cost = _objective(fresh, Int(gap_open), weights)
    if after_cost >= before_cost:
        return (before_cost, before_cost)
    rows = fresh^
    return (before_cost, after_cost)


def refine(
    rows: List[List[UInt8]],
    passes: Int,
    band: Int = 64,
    max_band: Int = 4096,
    anchored: Bool = True,
    anchor_costs: Costs = Costs.edit(),
    gap_open: Int = 0,
    seeded: Bool = False,
    group_band: Int = 0,
    period: Int = 0,
    weights: List[Int] = List[Int](),
) raises -> List[List[UInt8]]:
    """`rows` (one length, `-` gaps) refined on their sum of pairs for at most `passes` passes over the distinct
    sequences, each within `band` letters of its place, widened up to `max_band` while its best path touches the
    band's edge; with `gap_open`, every gap a pair opens costs that much more (see `_realign`). With `group_band`,
    each pass also realigns every split of a UPGMA tree over the rows, one side's profile against the other's
    within that many columns of their present merge (`_merge_groups`): moves of many rows at once, such as a
    stretch of copies shifted by a whole copy in several sequences, that no single row's move reaches. With a tandem
    repeat's `period`, the moves know its copies: a row's band reaches a copy and more either side, the group moves'
    two copies, and the group moves also split the rows by copy number (`_length_groups`). Identical
    sequences share their row, which weighs as many as spell it (each row counting its `weights`, 1 where empty): the
    objective is over every pair of the input rows. The rows keep their order and spell what they did."""
    var n = len(rows)
    if n < 2 or passes <= 0:
        return rows.copy()
    # One row per distinct sequence.
    var spelled = List[List[UInt8]]()
    var distinct = List[Int]()
    var row_of = List[Int](length=n, fill=-1)
    for r in range(n):
        var letters = List[UInt8]()
        for x in rows[r]:
            if x != GAP:
                letters.append(x)
        var found = -1
        for d in range(len(distinct)):
            if spelled[distinct[d]] == letters:
                found = d
                break
        spelled.append(letters^)
        if found < 0:
            row_of[r] = len(distinct)
            distinct.append(r)
        else:
            row_of[r] = found
    var work = List[List[UInt8]]()
    for d in distinct:
        work.append(rows[d].copy())
    if len(work) < 2:
        return rows.copy()
    var weight = List[Int](length=len(distinct), fill=0)
    for r in range(n):
        weight[row_of[r]] += _w(weights, r)
    # Columns gapped in every distinct row go first.
    var keep = List[Int]()
    for c in range(len(work[0])):
        for r in range(len(work)):
            if work[r][c] != GAP:
                keep.append(c)
                break
    for r in range(len(work)):
        var row = List[UInt8](capacity=len(keep))
        for c in keep:
            row.append(work[r][c])
        work[r] = row^
    var cost = List[Int32]()
    var moves = List[UInt8]()
    var centre = List[Int]()
    var spelled_distinct = List[List[UInt8]]()
    for d in distinct:
        spelled_distinct.append(spelled[d].copy())
    if not anchored:
        spelled_distinct.clear()
    var partners = Partners(spelled_distinct, anchor_costs, seeded)
    # A tandem repeat's copies: bands reaching a copy (a row) or two (a group) either side.
    var row_band = max(band, period + 2 * EDGE + 8) if period > 0 else band
    var groups_band = max(group_band, 2 * period + 16) if period > 0 else group_band
    var last = _objective(work, gap_open, weight)
    for _ in range(passes):
        for s in range(len(work)):
            centre.clear()
            if anchored:
                _anchored_centre(work, s, partners, centre)
            _ = _realign(work, s, row_band, max_band, centre, Int32(gap_open), cost, moves, List[Int](), False, weight)
        if group_band > 0 or period > 0:
            var groups = _tree_groups(work)
            if period > 0:
                groups.extend(_length_groups(work, period))
            for group in groups:
                _ = _merge_groups(
                    work,
                    group,
                    groups_band,
                    max(groups_band, min(max_band, 1024)),
                    Int32(gap_open),
                    cost,
                    moves,
                    weight,
                )
        var now = _objective(work, gap_open, weight)
        # A pass gaining under 0.1% ends the refinement: the later passes of a long run gain little.
        if now >= last or (last - now) * 1000 < last:
            break
        last = now
    var out = List[List[UInt8]]()
    for r in range(n):
        out.append(work[row_of[r]].copy())
    return out^


def place(
    mut rows: List[List[UInt8]],
    s: Int,
    centre: List[Int],
    forced: List[Int],
    band: Int = 64,
    max_band: Int = 4096,
    gap_open: Int = 1,
    weights: List[Int] = List[Int](),
) -> Tuple[Int, Int]:
    """Row `s` of `rows` realigned once, exactly, to the others' columns: the band centred on `centre` (sequence
    letters before each column, len(rows[0]) + 1 of them; empty for its present place), letter j held to column
    `forced[j]` where that is not -1 (`forced` empty for none), the other rows counting their `weights`. The new row
    is kept whatever it costs; returns `row_cost` before and after."""
    var cost = List[Int32]()
    var moves = List[UInt8]()
    return _realign(rows, s, band, max_band, centre, Int32(gap_open), cost, moves, forced, True, weights)
