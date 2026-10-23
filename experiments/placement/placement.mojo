"""Placing one sequence into an MSA through seed anchors: which seeds, combined how, place it closest to its floor.

usage: placement <msa.fa> <label> <out.tsv> <strategy> [<strategy> ...]
       placement --progressive <sequences.fa> <out.msa.fa> <strategy> [refine passes [group band[,max band]]]
       placement --refine <msa.fa> <out.msa.fa> <passes> [gap [anchored [band [group band]]]]

Progressive: the MSA grown a sequence at a time, longest first, each new sequence anchored by the strategy's seeds
against every sequence already in it and placed by the exact profile DP centred on its anchors (an MSA's first
sequence is its own row); then, with `refine passes`, plegma's refinement (-I passes,64,1). Identical sequences are
placed once and share their row.

Refine: an MSA's distinct rows through plegma's refinement (-I passes,band,gap with `anchored` 1 or 0, and the
group moves within `group band` columns when given), its objective before and after printed.

Leave one out: for every distinct row s of the MSA, its sequence is placed afresh into the other rows' columns, and
the placed row's exact cost against the others (`row_cost`: the columns where two rows differ, plus 1 for every gap
the pair opens) is written beside its present cost and its floor, the sum of its optimal pairwise costs to each
other row under the same costs (karan, mismatch 1, a gap of k letters 1 + k). A strategy is `seeds[;option ...]`:
- `seeds`: the seeding harness's levels (`mum51,mum31,mum21`, `smem31`, `ms21`, `rare8_21`, `combo8_21`, `all31`),
  run between the sequence and each other row's sequence through plegma's backbone; or `optimal` (every letter's
  partner in its optimal two-piece alignment to each row, what refine's anchors use now), `present` (no move) or
  `free` (realigned around its present place, no anchors).
- `stay`, `M=k:d` (repeat mask before chaining), `D=n` (the mask at the first n depths only), `T=n` (letters an
  anchor gives back at each end, 10): as in the seeding harness.
- `vote=row` (each row's seed one vote, the default) or `vote=len` (weighted by seed length); `share=F` (0.5): a
  letter is anchored where one column holds more than this share of the weight voted for it, by at least two rows;
  the anchors kept are the longest chain rising in letter and column.
- `hard`: the anchors fixed (each anchored letter held to its column), else the band centred on them;
  `band=N` (64).
The output: label, row, strategy, anchored letters, present cost, placed cost, floor, milliseconds.
"""

from std.sys import argv
from std.time import perf_counter_ns

from karan import Aligner, Costs, Mode, Ties
from sporos import read_fasta

from plegma.backbone import (
    ALL,
    Level,
    MUMS,
    MUMS_OR_SMEMS,
    RARE,
    RARE_OR_MAXIMAL,
    REFERENCE_UNIQUE,
    SMEMS,
    backbone,
)
from plegma.chain import ChainCosts
from plegma.graph import DELETION, INSERTION, cigar_runs
from plegma.refine import place, refine, row_cost

comptime GAP = UInt8(ord("-"))


struct Strategy(Copyable, Movable):
    var text: String
    var kind: String
    """"seeds", "optimal", "present" or "free"."""
    var levels: List[Level]
    var stay: Bool
    var trim: Int
    var mask_k: Int
    var mask_distance: Int
    var mask_depth: Int
    var by_length: Bool
    var share: Float64
    var hard: Bool
    var band: Int
    var sources: List[String]
    """The sources, best first: this one alone, or those its text joins with `>`."""

    def __init__(out self, text: String) raises:
        self.text = text
        self.kind = "seeds"
        self.levels = List[Level]()
        self.stay = False
        self.trim = 10
        self.mask_k = 0
        self.mask_distance = 0
        self.mask_depth = Int.MAX
        self.by_length = False
        self.share = 0.5
        self.hard = False
        self.band = 64
        self.sources = List[String]()
        if text.find(">") >= 0:
            # A cascade: each source its own strategy; the placement options (hard, band) from the first.
            for piece in text.split(">"):
                self.sources.append(String(piece.strip()))
            var first = Strategy(self.sources[0])
            self.kind = "cascade"
            self.hard = first.hard
            self.band = first.band
            return
        var parts = text.split(";")
        for k in range(1, len(parts)):
            var part = String(parts[k])
            if part == "stay":
                self.stay = True
            elif part == "hard":
                self.hard = True
            elif part.startswith("T="):
                self.trim = Int(String(part[byte=2:]))
            elif part.startswith("D="):
                self.mask_depth = Int(String(part[byte=2:]))
            elif part.startswith("M="):
                var halves = String(part[byte=2:]).split(":")
                self.mask_k = Int(String(halves[0]))
                self.mask_distance = Int(String(halves[1]))
            elif part == "vote=len":
                self.by_length = True
            elif part == "vote=row":
                self.by_length = False
            elif part.startswith("share="):
                self.share = Float64(String(part[byte=6:]))
            elif part.startswith("band="):
                self.band = Int(String(part[byte=5:]))
            else:
                raise Error(t"unknown option {part} in {text}")
        var head = String(parts[0])
        if head == "optimal" or head == "present" or head == "free":
            self.kind = head
            if head == "optimal":
                self.sources.append(text)
            return
        for item in head.split(","):
            var token = String(item)
            if token.startswith("mum"):
                self.levels.append(Level(MUMS, Int(String(token[byte=3:]))))
            elif token.startswith("smem"):
                self.levels.append(Level(SMEMS, Int(String(token[byte=4:]))))
            elif token.startswith("ms"):
                self.levels.append(Level(MUMS_OR_SMEMS, Int(String(token[byte=2:]))))
            elif token.startswith("ref"):
                self.levels.append(Level(REFERENCE_UNIQUE, Int(String(token[byte=3:]))))
            elif token.startswith("rare"):
                var halves = String(token[byte=4:]).split("_")
                self.levels.append(Level(RARE, Int(String(halves[1])), Int(String(halves[0]))))
            elif token.startswith("combo"):
                var halves = String(token[byte=5:]).split("_")
                self.levels.append(Level(RARE_OR_MAXIMAL, Int(String(halves[1])), Int(String(halves[0]))))
            elif token.startswith("all"):
                self.levels.append(Level(ALL, Int(String(token[byte=3:]))))
            else:
                raise Error(t"unknown seed level {token} in {text}")
        self.sources.append(text)


def _ungapped(row: List[UInt8]) -> List[UInt8]:
    var out = List[UInt8]()
    for x in row:
        if x != GAP:
            out.append(x)
    return out^


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


def _leading(
    st: Strategy,
    s: Int,
    q: List[UInt8],
    seen: List[List[UInt8]],
    column_of: List[List[Int]],
    mut aligner: Aligner,
    two_piece: Costs,
) raises -> List[Int]:
    """Each letter of the sequence's leading column by `st`'s votes from every other row (-1: not anchored)."""
    var n = len(q)
    var n_rows = len(seen)
    var vote_c = List[List[Int]]()
    var vote_w = List[List[Int]]()
    for _ in range(n):
        vote_c.append(List[Int]())
        vote_w.append(List[Int]())
    var voters = List[Int](length=n, fill=0)
    for r in range(n_rows):
        if r == s:
            continue
        if st.kind == "optimal":
            var found = aligner.align(Span(seen[r]), Span(q), two_piece, Mode.GLOBAL, ties=Ties.LEFT)
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
                        vote_c[j].append(column_of[r][i])
                        vote_w[j].append(1)
                        voters[j] += 1
                        i += 1
                        j += 1
        else:
            var anchors = backbone(
                Span(seen[r]),
                Span(q),
                st.levels,
                ChainCosts(substitution=0.1),
                0,
                st.stay,
                repeat_k=st.mask_k,
                repeat_distance=st.mask_distance,
                mask_depth=st.mask_depth,
            )
            for a in anchors:
                var lo = min(st.trim, a.length // 2)
                for t in range(lo, a.length - lo):
                    var j = a.query_pos + t
                    vote_c[j].append(column_of[r][a.ref_pos + t])
                    vote_w[j].append(a.length if st.by_length else 1)
                    voters[j] += 1
    var leading = List[Int](length=n, fill=-1)
    for j in range(n):
        if voters[j] < 2:
            continue
        var total = 0
        var best_c = -1
        var best_w = 0
        for x in range(len(vote_c[j])):
            total += vote_w[j][x]
            var c = vote_c[j][x]
            var w = 0
            for y in range(len(vote_c[j])):
                if vote_c[j][y] == c:
                    w += vote_w[j][y]
            if w > best_w:
                best_w = w
                best_c = c
        if Float64(best_w) > st.share * Float64(total):
            leading[j] = best_c
    return leading^


def _heaviest_chain(leading: List[Int], weight: List[Int], columns: Int) -> List[Int]:
    """The letters whose leading columns form the heaviest chain rising in letter and column (a Fenwick tree of the
    best chain ending before each column)."""
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
        # The best chain over columns < c.
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


def _progressive(args: List[String]) raises:
    var records = read_fasta(args[2])
    var st = Strategy(args[4])
    var parts = List[Strategy]()
    for text in st.sources:
        parts.append(Strategy(text))
    var passes = Int(args[5]) if len(args) > 5 else 0
    # The refinement's group band and its max band: `group[,max band]` (0 and 512 by default).
    var group_band = 0
    var max_band = 512
    if len(args) > 6:
        var parts_text = args[6].split(",")
        group_band = Int(parts_text[0])
        if len(parts_text) > 1:
            max_band = Int(parts_text[1])
    var seqs = List[List[UInt8]]()
    var names = List[String]()
    for r in range(len(records)):
        var letters = List[UInt8](records.seq(r))
        for k in range(len(letters)):
            if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
                letters[k] -= 32
        seqs.append(letters^)
        names.append(records.names[r])
    # Distinct sequences, longest first.
    var distinct = List[Int]()
    var of = List[Int](length=len(seqs), fill=-1)
    for r in range(len(seqs)):
        for d in range(len(distinct)):
            if seqs[distinct[d]] == seqs[r]:
                of[r] = d
                break
        if of[r] < 0:
            of[r] = len(distinct)
            distinct.append(r)
    var keys = List[UInt64]()
    for d in range(len(distinct)):
        keys.append((UInt64((1 << 32) - len(seqs[distinct[d]])) << 32) | UInt64(d))
    sort(keys)
    var order = List[Int]()
    for key in keys:
        order.append(Int(key & 0xFFFFFFFF))
    var aligner = Aligner()
    var two_piece = Costs.two_piece(4, 4, 2, 24, 1)
    # rows[k]: the MSA row of distinct sequence order[k]; placed: their sequences.
    var rows = List[List[UInt8]]()
    var placed = List[List[UInt8]]()
    for step in range(len(order)):
        var q = seqs[distinct[order[step]]].copy()
        if step == 0:
            rows.append(q.copy())
            placed.append(q^)
            continue
        var columns = len(rows[0])
        var n = len(q)
        # The new row first takes columns of its own after the others' (every letter inserted).
        for r in range(len(rows)):
            for _ in range(n):
                rows[r].append(GAP)
        var mine = List[UInt8](length=columns, fill=GAP)
        mine.extend(q.copy())
        rows.append(mine^)
        placed.append(q.copy())
        var s = len(rows) - 1
        var column_of = List[List[Int]]()
        for r in range(len(rows)):
            var cols = List[Int]()
            for c in range(len(rows[r])):
                if rows[r][c] != GAP:
                    cols.append(c)
            column_of.append(cols^)
        var leading = List[Int](length=n, fill=-1)
        var weight = List[Int](length=n, fill=0)
        if st.kind != "free" and st.kind != "present":
            var rank_weight = 1
            for _ in range(len(parts)):
                rank_weight *= 1000
            for part in parts:
                rank_weight //= 1000
                var got = _leading(part, s, q, placed, column_of, aligner, two_piece)
                for j in range(n):
                    if leading[j] < 0 and got[j] >= 0:
                        leading[j] = got[j]
                        weight[j] = rank_weight
            # A single placed partner votes alone: accept its anchors too.
            if step == 1:
                var got = _leading_one(parts[0], q, placed[0], column_of[0], aligner, two_piece)
                for j in range(n):
                    if leading[j] < 0 and got[j] >= 0:
                        leading[j] = got[j]
                        weight[j] = 1
        var total_columns = len(rows[0])
        var chain_j = _heaviest_chain(leading, weight, total_columns)
        var pc = List[Int]()
        var pj = List[Int]()
        for j in chain_j:
            pc.append(leading[j])
            pj.append(j)
        var centre = _centre_of(pc, pj, total_columns, n)
        if len(chain_j) == 0:
            # No anchor: the diagonal through the others' columns.
            for c in range(total_columns + 1):
                centre[c] = min(n, c * n // max(1, columns))
        _ = place(rows, s, centre, List[Int](), st.band, 4096, 1)
    var final_rows = List[List[UInt8]]()
    for k in range(len(order)):
        final_rows.append(rows[k].copy())
    if passes > 0:
        final_rows = refine(final_rows, passes, 64, max_band, True, two_piece, 1, False, group_band)
    var row_of = List[Int](length=len(distinct), fill=-1)
    for k in range(len(order)):
        row_of[order[k]] = k
    with open(args[3], "w") as f:
        for r in range(len(seqs)):
            f.write(">", names[r], "\n", String(unsafe_from_utf8=final_rows[row_of[of[r]]].copy()), "\n")


def _leading_one(
    st: Strategy, q: List[UInt8], partner: List[UInt8], column_of: List[Int], mut aligner: Aligner, two_piece: Costs
) raises -> List[Int]:
    """The columns one partner's seeds (or optimal alignment) give the sequence's letters, without a vote."""
    var leading = List[Int](length=len(q), fill=-1)
    if st.kind == "optimal":
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
                    leading[j] = column_of[i]
                    i += 1
                    j += 1
        return leading^
    var anchors = backbone(
        Span(partner),
        Span(q),
        st.levels,
        ChainCosts(substitution=0.1),
        0,
        st.stay,
        repeat_k=st.mask_k,
        repeat_distance=st.mask_distance,
        mask_depth=st.mask_depth,
    )
    for a in anchors:
        var lo = min(st.trim, a.length // 2)
        for t in range(lo, a.length - lo):
            leading[a.query_pos + t] = column_of[a.ref_pos + t]
    return leading^


def _refine_msa(args: List[String]) raises:
    """An MSA's distinct rows through plegma's refinement, every input name written with its row."""
    var records = read_fasta(args[2])
    var passes = Int(args[4])
    var gap = Int(args[5]) if len(args) > 5 else 1
    var anchored = Int(args[6]) != 0 if len(args) > 6 else True
    var band = Int(args[7]) if len(args) > 7 else 64
    var group_band = Int(args[8]) if len(args) > 8 else 0
    var rows = List[List[UInt8]]()
    var seen = List[List[UInt8]]()
    var of = List[Int]()
    var names = List[String]()
    for r in range(len(records)):
        var row = List[UInt8](records.seq(r))
        for k in range(len(row)):
            if row[k] >= UInt8(ord("a")) and row[k] <= UInt8(ord("z")):
                row[k] -= 32
            if row[k] == UInt8(ord(".")):
                row[k] = GAP
        names.append(records.names[r])
        var spelled = _ungapped(row)
        var at = -1
        for d in range(len(seen)):
            if seen[d] == spelled:
                at = d
                break
        if at < 0:
            at = len(seen)
            seen.append(spelled^)
            rows.append(row^)
        of.append(at)
    var keep = List[Int]()
    for c in range(len(rows[0])):
        for r in range(len(rows)):
            if rows[r][c] != GAP:
                keep.append(c)
                break
    for r in range(len(rows)):
        var row = List[UInt8](capacity=len(keep))
        for c in keep:
            row.append(rows[r][c])
        rows[r] = row^
    var t0 = perf_counter_ns()
    var before = 0
    for s in range(len(rows)):
        before += row_cost(rows, s, gap)
    if passes > 0:
        rows = refine(rows, passes, band, 4096, anchored, Costs.two_piece(4, 4, 2, 24, 1), gap, False, group_band)
    var after = 0
    for s in range(len(rows)):
        after += row_cost(rows, s, gap)
    print("rows", len(rows), "objective", before // 2, "->", after // 2, "ms", (perf_counter_ns() - t0) // 1_000_000)
    with open(args[3], "w") as f:
        for r in range(len(names)):
            f.write(">", names[r], "\n", String(unsafe_from_utf8=rows[of[r]].copy()), "\n")


def main() raises:
    var raw = argv()
    if len(raw) > 1 and (String(raw[1]) == "--progressive" or String(raw[1]) == "--refine"):
        var list = List[String]()
        for k in range(len(raw)):
            list.append(String(raw[k]))
        if list[1] == "--progressive":
            _progressive(list)
        else:
            _refine_msa(list)
        return
    var args = argv()
    if len(args) < 5:
        print("usage: placement <msa.fa> <label> <out.tsv> <strategy> [<strategy> ...]")
        return
    var records = read_fasta(String(args[1]))
    var label = String(args[2])
    var strategies = List[Strategy]()
    var sources_of = List[List[Strategy]]()
    for k in range(4, len(args)):
        var st = Strategy(String(args[k]))
        var parts = List[Strategy]()
        for text in st.sources:
            parts.append(Strategy(text))
        strategies.append(st^)
        sources_of.append(parts^)
    # Distinct rows, columns gapped in every one dropped.
    var rows = List[List[UInt8]]()
    var names = List[String]()
    var seen = List[List[UInt8]]()
    for r in range(len(records)):
        var row = List[UInt8](records.seq(r))
        for k in range(len(row)):
            if row[k] >= UInt8(ord("a")) and row[k] <= UInt8(ord("z")):
                row[k] -= 32
            if row[k] == UInt8(ord(".")):
                row[k] = GAP
        var spelled = _ungapped(row)
        var dup = False
        for x in seen:
            if x == spelled:
                dup = True
                break
        if dup:
            continue
        seen.append(spelled^)
        rows.append(row^)
        names.append(records.names[r])
    var keep = List[Int]()
    for c in range(len(rows[0])):
        for r in range(len(rows)):
            if rows[r][c] != GAP:
                keep.append(c)
                break
    for r in range(len(rows)):
        var row = List[UInt8](capacity=len(keep))
        for c in keep:
            row.append(rows[r][c])
        rows[r] = row^
    var n_rows = len(rows)
    var columns = len(rows[0])
    var aligner = Aligner()
    var floor_costs = Costs.affine(1, 1, 1)
    var two_piece = Costs.two_piece(4, 4, 2, 24, 1)
    var out = String()
    # Each row's letters' columns.
    var column_of = List[List[Int]]()
    for r in range(n_rows):
        var cols = List[Int]()
        for c in range(columns):
            if rows[r][c] != GAP:
                cols.append(c)
        column_of.append(cols^)
    for s in range(n_rows):
        var q = seen[s].copy()
        var n = len(q)
        var present = row_cost(rows, s, 1)
        var floor = 0
        for r in range(n_rows):
            if r != s:
                floor += aligner.distance(Span(q), Span(seen[r]), floor_costs)
        for x in range(len(strategies)):
            ref st = strategies[x]
            var t0 = perf_counter_ns()
            var work = rows.copy()
            var anchored = 0
            var result = (present, present)
            if st.kind == "present":
                pass
            elif st.kind == "free":
                result = place(work, s, List[Int](), List[Int](), st.band, 4096, 1)
            else:
                # Each source's leading columns; a letter takes the first source's that anchors it, and the chain
                # is the heaviest rising one, each source outweighing every later one.
                var leading = List[Int](length=n, fill=-1)
                var weight = List[Int](length=n, fill=0)
                ref parts = sources_of[x]
                var rank_weight = 1
                for _ in range(len(parts)):
                    rank_weight *= 1000
                for part in parts:
                    rank_weight //= 1000
                    var mine = _leading(part, s, q, seen, column_of, aligner, two_piece)
                    for j in range(n):
                        if leading[j] < 0 and mine[j] >= 0:
                            leading[j] = mine[j]
                            weight[j] = rank_weight
                var chain_j = _heaviest_chain(leading, weight, columns)
                var pc = List[Int]()
                var pj = List[Int]()
                var forced = List[Int]()
                if st.hard:
                    forced = List[Int](length=n, fill=-1)
                for j in chain_j:
                    pc.append(leading[j])
                    pj.append(j)
                    if st.hard:
                        forced[j] = leading[j]
                anchored = len(chain_j)
                var centre = _centre_of(pc, pj, columns, n)
                result = place(work, s, centre, forced, st.band, 4096, 1)
            var ms = Float64(perf_counter_ns() - t0) / 1e6
            out.write(
                label,
                "\t",
                names[s],
                "\t",
                st.text,
                "\t",
                anchored,
                "\t",
                result[0],
                "\t",
                result[1],
                "\t",
                floor,
                "\t",
                ms,
                "\n",
            )
    with open(String(args[3]), "w") as f:
        f.write("label\trow\tstrategy\tanchored\tpresent\tplaced\tfloor\tms\n")
        f.write(out)
