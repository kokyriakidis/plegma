"""Seed anchors against the truth: which anchors pair letters no optimal alignment pairs, and what marks them.

usage: anchors <sequences.fa> <period> <max pairs> <out.tsv> <plan> [<plan> ...]

For an even sample of at most `max pairs` pairs of distinct sequences (the longer as the partner, as placement
seeds a new sequence against a longer placed one), each plan seeds the pair as plegma's backbone does, traced: every
anchor with the depth of the recursion and the level that found it, its length, its string's copies in each piece,
the window it was found in, whether the partner's repeat mask covers it, how many anchors its window chained, and
its diagonal against the window's two ends. Each anchored letter pair (past the trim) is checked against every
optimal unit-cost alignment of the pair, vntr-eval's metric: on one exactly when its forward cost, its own and its
backward cost sum to the optimum. A letter pair on none is off; its offset is the distance to the nearest partner
an optimal alignment gives the letter, in whole periods (`copies`), within two letters (`near`) or else (`other`).

A plan is `levels[;option ...]`: levels `mum51,mum31,mum21` (MUMs), `smem31` (SMEMs), `ref31` (unique in the
partner); options `stay` (a level again inside its own gaps until one yields nothing), `M=k:d` (no anchor in the
partner's tandem repeats: k-mers recurring within d letters), `D=n` (that mask at the first n depths only), `T=n`
(letters an anchor gives back at each end, 10), `B=n` (no seeding inside a window of n letters or fewer: left to the
exact aligner), `ends=e` (past the mask's depths, an anchor in the repeat is kept only if its diagonal lies within e
letters of one of its window's ends: the window's length difference all after it or all before it, one indel).
Or `optimal` (the optimal two-piece alignment, 4,4,2,24,1, what placement's last source uses) or `edit` (the
optimal unit-cost one): every matched letter pair, as runs of matches. Or `copies[;min=L]`: the flanks' MUMs (51,
31, 21, masked in the repeats at every depth), then inside the arrays (`units.repeat_at` at the period, cut into
copies against the region's unit, the longest sequence's consensus) the copies paired by `units.pair_copies`, and
inside each pair of copies the runs of at least L (1) matching letters of their optimal unit-cost alignment.

Being off an optimal alignment is a strict test in a repeat: an alignment a fraction of a percent worse may pair
most letters a copy apart. So each plan is also charged what its anchors force: the pair's best alignment through
every anchored letter pair (each piece between consecutive anchors aligned optimally, as vntr-eval scores an MSA's
shared letters) minus its optimum, summed over the pairs, against the summed optima.

Writes one TSV line per anchor (plan, pair, its features, letters, off letters by offset class) and prints per plan
the anchors, the letters and how many are off, the excess the anchors force, and per depth the anchors and letters
off.
"""

from std.sys import argv

from karan import Aligner, Costs, Mode, Ties
from sporos import read_fasta

from plegma.anchors import Anchor, find_anchors
from plegma.backbone import Level, MUMS, REFERENCE_UNIQUE, SMEMS, RepeatMask
from plegma.chain import ChainCosts, chain, without_overlaps
from plegma.graph import DELETION, INSERTION, cigar_runs
from plegma.units import Copies, consensus_unit, copies, pair_copies, repeat_at


struct Plan(Copyable, Movable):
    var text: String
    var levels: List[Level]
    var stay: Bool
    var trim: Int
    var mask_k: Int
    var mask_distance: Int
    var mask_depth: Int
    var budget: Int
    var ends: Int
    var min_run: Int
    var kind: String
    """"seeds", "optimal", "edit" or "copies"."""

    def __init__(out self, text: String) raises:
        self.text = text
        self.ends = -1
        self.min_run = 1
        self.kind = "seeds"
        self.levels = List[Level]()
        self.stay = False
        self.trim = 10
        self.mask_k = 0
        self.mask_distance = 0
        self.mask_depth = Int.MAX
        self.budget = 0
        var parts = text.split(";")
        for k in range(1, len(parts)):
            var part = String(parts[k])
            if part == "stay":
                self.stay = True
            elif part.startswith("T="):
                self.trim = Int(String(part[byte=2:]))
            elif part.startswith("D="):
                self.mask_depth = Int(String(part[byte=2:]))
            elif part.startswith("B="):
                self.budget = Int(String(part[byte=2:]))
            elif part.startswith("ends="):
                self.ends = Int(String(part[byte=5:]))
            elif part.startswith("min="):
                self.min_run = Int(String(part[byte=4:]))
            elif part.startswith("M="):
                var halves = String(part[byte=2:]).split(":")
                self.mask_k = Int(String(halves[0]))
                self.mask_distance = Int(String(halves[1]))
            else:
                raise Error(t"unknown option {part} in {text}")
        if String(parts[0]) == "optimal" or String(parts[0]) == "edit":
            self.kind = String(parts[0])
            self.levels.append(Level(MUMS, 1))
            return
        if String(parts[0]) == "copies":
            # The flanks' MUMs, masked at every depth; the copies' anchors are added to them.
            self.kind = "copies"
            self.trim = 0
            self.mask_k = 12
            self.mask_distance = 500
            for k in [51, 31, 21]:
                self.levels.append(Level(MUMS, k))
            self.stay = True
            return
        for item in String(parts[0]).split(","):
            var token = String(item)
            if token.startswith("mum"):
                self.levels.append(Level(MUMS, Int(String(token[byte=3:]))))
            elif token.startswith("smem"):
                self.levels.append(Level(SMEMS, Int(String(token[byte=4:]))))
            elif token.startswith("ref"):
                self.levels.append(Level(REFERENCE_UNIQUE, Int(String(token[byte=3:]))))
            else:
                raise Error(t"unknown level {token} in {text}")


@fieldwise_init
struct Traced(Copyable, Movable):
    var ref_pos: Int
    var query_pos: Int
    var length: Int
    var depth: Int
    var level: Int
    var ref_copies: Int
    var query_copies: Int
    var win_ref: Int
    var win_query: Int
    var in_repeat: Bool
    var siblings: Int
    var diag_start: Int
    """The anchor's diagonal against its window's start: query letters before it minus partner letters."""
    var diag_end: Int
    """Likewise against the window's end: query letters after it minus partner letters."""


def _trace(
    a: List[UInt8],
    b: List[UInt8],
    r0: Int,
    r1: Int,
    q0: Int,
    q1: Int,
    plan: Plan,
    level: Int,
    depth: Int,
    repeats: RepeatMask,
    mut out: List[Traced],
) raises:
    """`backbone`'s recursion over window [r0, r1) x [q0, q1), every anchor kept traced."""
    if level >= len(plan.levels) or max(r1 - r0, q1 - q0) <= plan.budget:
        return
    var current = plan.levels[level]
    var found = current.select(find_anchors(Span(a)[r0:r1], Span(b)[q0:q1], current.min_len, current.counted_to()))
    if repeats.k > 0 and depth < plan.mask_depth:
        var outside = List[Anchor]()
        for x in found:
            if not repeats.covers(r0 + x.ref_pos, r0 + x.ref_end()):
                outside.append(x)
        found = outside^
    elif repeats.k > 0 and plan.ends >= 0:
        # Inside the repeat, only anchors leaving the window's length difference to one side of them.
        var kept = List[Anchor]()
        for x in found:
            var start = x.query_pos - x.ref_pos
            var end = ((q1 - q0) - x.query_end()) - ((r1 - r0) - x.ref_end())
            if not repeats.covers(r0 + x.ref_pos, r0 + x.ref_end()) or abs(start) <= plan.ends or abs(end) <= plan.ends:
                kept.append(x)
        found = kept^
    var picked = List[Anchor]()
    if len(found) <= 50_000:
        picked = without_overlaps(chain(found, r1 - r0, q1 - q0, ChainCosts(substitution=0.1), 0).anchors)
    var next = level if plan.stay else level + 1
    if len(picked) == 0:
        _trace(a, b, r0, r1, q0, q1, plan, level + 1, depth, repeats, out)
        return
    var ref_at = 0
    var query_at = 0
    for x in picked:
        _trace(a, b, r0 + ref_at, r0 + x.ref_pos, q0 + query_at, q0 + x.query_pos, plan, next, depth + 1, repeats, out)
        out.append(
            Traced(
                r0 + x.ref_pos,
                q0 + x.query_pos,
                x.length,
                depth,
                level,
                x.ref_copies,
                x.query_copies,
                r1 - r0,
                q1 - q0,
                repeats.covers(r0 + x.ref_pos, r0 + x.ref_end()),
                len(picked),
                x.query_pos - x.ref_pos,
                ((q1 - q0) - x.query_end()) - ((r1 - r0) - x.ref_end()),
            )
        )
        ref_at = x.ref_end()
        query_at = x.query_end()
    _trace(a, b, r0 + ref_at, r1, q0 + query_at, q1, plan, next, depth + 1, repeats, out)


def _aligned(a: List[UInt8], b: List[UInt8], unit: Bool, mut out: List[Traced]) raises:
    """The matched letter pairs of the pair's optimal alignment (two-piece, or unit costs), as runs of matches."""
    var aligner = Aligner()
    var costs = Costs.edit() if unit else Costs.two_piece(4, 4, 2, 24, 1)
    var found = aligner.align(Span(a), Span(b), costs, Mode.GLOBAL, ties=Ties.LEFT)
    var runs = cigar_runs(found.cigar)
    var i = 0
    var j = 0
    var run_i = -1
    var run_j = -1
    var run = 0

    @parameter
    def close():
        if run > 0:
            out.append(Traced(run_i, run_j, run, 0, 0, 1, 1, len(a), len(b), False, 0, run_j - run_i, 0))

    for x in range(len(runs.ops)):
        for _ in range(runs.counts[x]):
            var op = runs.ops[x]
            if op == DELETION:
                close()
                run = 0
                i += 1
            elif op == INSERTION:
                close()
                run = 0
                j += 1
            else:
                if a[i] == b[j]:
                    if run == 0:
                        run_i = i
                        run_j = j
                    run += 1
                else:
                    close()
                    run = 0
                i += 1
                j += 1
    close()


def _copy_anchors(
    a: List[UInt8], b: List[UInt8], period: Int, unit: List[UInt8], min_run: Int, mut out: List[Traced]
) raises:
    """Inside the two arrays, the copies paired (`pair_copies`) and the matching runs of each pair's optimal
    unit-cost alignment, at least `min_run` long; those crossing an anchor already in `out` (the flanks') dropped."""
    var ra = repeat_at(a, period)
    var rb = repeat_at(b, period)
    if ra.period == 0 or rb.period == 0:
        return
    var ca = copies(a, ra, unit)
    var cb = copies(b, rb, unit)
    var aligner = Aligner()
    var pairs = pair_copies(a, ca, b, cb, aligner)
    var flanks = len(out)
    for pr in pairs:
        var i0 = ca.starts[pr[0]]
        var j0 = cb.starts[pr[1]]
        var found = aligner.align(
            Span(a)[i0 : ca.ends[pr[0]]], Span(b)[j0 : cb.ends[pr[1]]], Costs.edit(), Mode.GLOBAL, ties=Ties.LEFT
        )
        var runs = cigar_runs(found.cigar)
        var i = i0
        var j = j0
        var run_i = -1
        var run_j = -1
        var run = 0
        for x in range(len(runs.ops) + 1):
            var count = runs.counts[x] if x < len(runs.ops) else 1
            for _ in range(count):
                var matched = (
                    x < len(runs.ops) and runs.ops[x] != DELETION and runs.ops[x] != INSERTION and a[i] == b[j]
                )
                if matched:
                    if run == 0:
                        run_i = i
                        run_j = j
                    run += 1
                elif run >= min_run:
                    # Kept unless it crosses a flank anchor (not before or after it on both sequences).
                    var crosses = False
                    for f in range(flanks):
                        ref fa = out[f]
                        var before = run_i + run <= fa.ref_pos and run_j + run <= fa.query_pos
                        var after = run_i >= fa.ref_pos + fa.length and run_j >= fa.query_pos + fa.length
                        if not before and not after:
                            crosses = True
                            break
                    if not crosses:
                        out.append(
                            Traced(run_i, run_j, run, 1, 0, 1, 1, len(a), len(b), True, len(pairs), run_j - run_i, 0)
                        )
                    run = 0
                else:
                    run = 0
                if x >= len(runs.ops):
                    break
                var op = runs.ops[x]
                if op == DELETION:
                    i += 1
                elif op == INSERTION:
                    j += 1
                else:
                    i += 1
                    j += 1


def main() raises:
    var args = argv()
    if len(args) < 6:
        print("usage: anchors <sequences.fa> <period> <max pairs> <out.tsv> <plan> [<plan> ...]")
        return
    var records = read_fasta(String(args[1]))
    var period = Int(String(args[2]))
    var max_pairs = Int(String(args[3]))
    var plans = List[Plan]()
    for k in range(5, len(args)):
        plans.append(Plan(String(args[k])))
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        var dup = False
        for x in seqs:
            if x == letters:
                dup = True
                break
        if not dup:
            seqs.append(letters^)
    # The region's unit: the longest sequence's array at the period, its consensus.
    var unit = List[UInt8]()
    if period > 0:
        var longest = 0
        for k in range(1, len(seqs)):
            if len(seqs[k]) > len(seqs[longest]):
                longest = k
        var top = repeat_at(seqs[longest], period)
        if top.period > 0:
            unit = consensus_unit(seqs[longest], top.start, top.end, period)
    var pairs = List[Tuple[Int, Int]]()
    for x in range(len(seqs)):
        for y in range(x + 1, len(seqs)):
            pairs.append((x, y) if len(seqs[x]) >= len(seqs[y]) else (y, x))
    var chosen = List[Tuple[Int, Int]]()
    var step = max(1, len(pairs) // max(1, max_pairs))
    var t = 0
    while t < len(pairs) and len(chosen) < max_pairs:
        chosen.append(pairs[t])
        t += step
    var n_plans = len(plans)
    comptime DEPTHS = 8
    # Per plan and depth (the last one gathering deeper ones): anchors, off anchors, letters, off letters.
    var tally = List[Int](length=n_plans * DEPTHS * 4, fill=0)
    var by_class = List[Int](length=n_plans * 3, fill=0)
    var forced = List[Int](length=n_plans, fill=0)
    var optimum = 0
    var aligner = Aligner()
    var text = String()
    text.write(
        "plan\ta\tb\tlen_a\tlen_b\tdepth\tlevel\tk\tlength\tref_copies\tquery_copies\twin_ref\twin_query\tin_repeat",
        "\tsiblings\tdiag_start\tdiag_end\tletters\toff\toff_near\toff_copies\toff_other\n",
    )
    var done = 0
    for pq in chosen:
        ref a = seqs[pq[0]]
        ref b = seqs[pq[1]]
        var n = len(a)
        var m = len(b)
        if (n + 1) * (m + 1) > 150_000_000:
            continue
        done += 1
        var w = m + 1
        var fwd = List[Int32](length=(n + 1) * w, fill=0)
        for j in range(w):
            fwd[j] = Int32(j)
        for i in range(1, n + 1):
            fwd[i * w] = Int32(i)
            var x = a[i - 1]
            for j in range(1, w):
                var d = fwd[(i - 1) * w + j - 1] + (Int32(0) if x == b[j - 1] else Int32(1))
                fwd[i * w + j] = min(d, min(fwd[(i - 1) * w + j] + 1, fwd[i * w + j - 1] + 1))
        var best = fwd[n * w + m]
        optimum += Int(best)
        for pi in range(n_plans):
            ref plan = plans[pi]
            var repeats = RepeatMask(Span(a), plan.mask_k, plan.mask_distance)
            var traced = List[Traced]()
            if plan.kind == "seeds":
                _trace(a, b, 0, n, 0, m, plan, 0, 0, repeats, traced)
            elif plan.kind == "copies":
                _trace(a, b, 0, n, 0, m, plan, 0, 0, repeats, traced)
                if period > 0:
                    _copy_anchors(a, b, period, unit, plan.min_run, traced)
            else:
                _aligned(a, b, plan.kind == "edit", traced)
            # Each anchored letter pair, by partner letter: (query letter, anchor).
            var at = List[List[Int]]()
            var owner = List[List[Int]]()
            for _ in range(n):
                at.append(List[Int]())
                owner.append(List[Int]())
            var letters = List[Int](length=len(traced), fill=0)
            var off = List[Int](length=len(traced), fill=0)
            var near = List[Int](length=len(traced), fill=0)
            var copies = List[Int](length=len(traced), fill=0)
            var other = List[Int](length=len(traced), fill=0)
            for k in range(len(traced)):
                ref an = traced[k]
                var lo = min(plan.trim, an.length // 2)
                for q in range(lo, an.length - lo):
                    at[an.ref_pos + q].append(an.query_pos + q)
                    owner[an.ref_pos + q].append(k)
            var next_row = List[Int32](length=w, fill=0)
            var row = List[Int32](length=w, fill=0)
            for j in range(w):
                next_row[j] = Int32(m - j)
            for i in reversed(range(n)):
                row[m] = Int32(n - i)
                var x = a[i]
                for j in reversed(range(m)):
                    var d = next_row[j + 1] + (Int32(0) if x == b[j] else Int32(1))
                    row[j] = min(d, min(next_row[j] + 1, row[j + 1] + 1))
                for e in range(len(at[i])):
                    var j = at[i][e]
                    var k = owner[i][e]
                    letters[k] += 1
                    if fwd[i * w + j] + (Int32(0) if x == b[j] else Int32(1)) + next_row[j + 1] == best:
                        continue
                    off[k] += 1
                    var offset = m
                    for dj in range(1, m):
                        var hit = False
                        for jj in [j - dj, j + dj]:
                            if jj >= 0 and jj < m:
                                if fwd[i * w + jj] + (Int32(0) if x == b[jj] else Int32(1)) + next_row[jj + 1] == best:
                                    hit = True
                        if hit:
                            offset = dj
                            break
                    if offset <= 2:
                        near[k] += 1
                    elif period > 0 and abs(offset - ((offset + period // 2) // period) * period) <= 2:
                        copies[k] += 1
                    else:
                        other[k] += 1
                next_row = row.copy()
            # What the anchors force: the pieces between consecutive anchored letters (kept while they rise on both
            # sequences, by partner letter) aligned optimally.
            var keys = List[UInt64]()
            for k in range(len(traced)):
                keys.append((UInt64(traced[k].ref_pos) << 32) | UInt64(k))
            sort(keys)
            var cost = 0
            var pi_at = 0
            var pj_at = 0
            for key in keys:
                ref an = traced[Int(key & 0xFFFFFFFF)]
                var lo = min(plan.trim, an.length // 2)
                var i0 = an.ref_pos + lo
                var j0 = an.query_pos + lo
                var len_kept = an.length - 2 * lo
                if len_kept <= 0 or i0 < pi_at or j0 < pj_at:
                    continue
                cost += aligner.distance(Span(a)[pi_at:i0], Span(b)[pj_at:j0])
                pi_at = i0 + len_kept
                pj_at = j0 + len_kept
            cost += aligner.distance(Span(a)[pi_at:n], Span(b)[pj_at:m])
            forced[pi] += cost - Int(best)
            for k in range(len(traced)):
                ref an = traced[k]
                var d = min(an.depth, DEPTHS - 1)
                var base = (pi * DEPTHS + d) * 4
                tally[base] += 1
                tally[base + 1] += 1 if off[k] > 0 else 0
                tally[base + 2] += letters[k]
                tally[base + 3] += off[k]
                by_class[pi * 3] += near[k]
                by_class[pi * 3 + 1] += copies[k]
                by_class[pi * 3 + 2] += other[k]
                text.write(
                    plan.text,
                    "\t",
                    pq[0],
                    "\t",
                    pq[1],
                    "\t",
                    n,
                    "\t",
                    m,
                    "\t",
                    an.depth,
                    "\t",
                    an.level,
                    "\t",
                    plan.levels[an.level].min_len,
                    "\t",
                    an.length,
                    "\t",
                    an.ref_copies,
                    "\t",
                    an.query_copies,
                    "\t",
                    an.win_ref,
                    "\t",
                    an.win_query,
                    "\t",
                    1 if an.in_repeat else 0,
                    "\t",
                    an.siblings,
                    "\t",
                    an.diag_start,
                    "\t",
                    an.diag_end,
                    "\t",
                    letters[k],
                    "\t",
                    off[k],
                    "\t",
                    near[k],
                    "\t",
                    copies[k],
                    "\t",
                    other[k],
                    "\n",
                )
    with open(String(args[4]), "w") as f:
        f.write(text)
    print("pairs", done, "of", len(pairs), "period", period)
    for pi in range(n_plans):
        var anchors = 0
        var off_anchors = 0
        var letters = 0
        var off_letters = 0
        for d in range(DEPTHS):
            var base = (pi * DEPTHS + d) * 4
            anchors += tally[base]
            off_anchors += tally[base + 1]
            letters += tally[base + 2]
            off_letters += tally[base + 3]
        print(
            plans[pi].text,
            "| anchors",
            anchors,
            "off",
            off_anchors,
            "| letters",
            letters,
            "off",
            off_letters,
            "| off by: near",
            by_class[pi * 3],
            "copies",
            by_class[pi * 3 + 1],
            "other",
            by_class[pi * 3 + 2],
            "| forced excess",
            forced[pi],
            "of optima",
            optimum,
        )
        for d in range(DEPTHS):
            var base = (pi * DEPTHS + d) * 4
            if tally[base] > 0:
                print(
                    "    depth",
                    d,
                    "anchors",
                    tally[base],
                    "off",
                    tally[base + 1],
                    "letters",
                    tally[base + 2],
                    "off",
                    tally[base + 3],
                )
