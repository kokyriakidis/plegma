"""Whether an MSA aligns a tandem repeat's copies as an optimal alignment does: copy by copy, pair by pair.

usage: structure <sequences.fa> <msa.fa> <period> <max pairs> [<msa.fa> ...] [--affine]
       structure <sequences.fa> <msa.fa> <period> --pair <a> <b> [--affine]

With `--affine`, the truth is the pairs' optimal gap-affine alignments (vntr-eval's affine costs: 4 per mismatch, a
gap of k letters 6 + 2k), which take a copy-number difference as few long gaps, not the unit-cost ones.

Every distinct sequence's array (`units.repeat_at` at the period) is cut into copies against the region's unit (the
longest sequence's consensus, `units.copies`). For an even sample of at most `max pairs` pairs (the longer as a, as
the eval's pairs are symmetric), the pair's forward and backward unit-cost DPs give each letter of a the letters of b
an optimal alignment pairs it with; a copy of a may pair with the copies of b holding those (its optimal range), or
be gapped whole where no optimal alignment pairs any of its letters. The MSA pairs a copy with the copy of b most of
its letters share columns with (none: under half of them face letters of b). A copy is
- `right`: the MSA's partner in its optimal range, or gapped where an optimal alignment gaps it;
- `off k`: its partner k copies outside the range (a copy shifted by k);
- `gapped`: the MSA gaps it where every optimal alignment pairs it; `paired`: the MSA pairs it where they gap it.
Tallied by where the copy sits in a's array (first tenth, middle, last tenth), by whether its string is a rare
variant (in two sequences or fewer, the copies carrying homology of their own) or common, and the flanks' letters
(outside the arrays) on or off an optimal alignment; then per sequence, the share of its copies right. Each MSA is
audited on the same pairs.
"""

from std.sys import argv

from karan import Aligner
from sporos import read_fasta

from plegma.units import Copies, consensus_unit, copies, repeat_at

comptime GAP = UInt8(ord("-"))


def _upper(var letters: List[UInt8]) -> List[UInt8]:
    for k in range(len(letters)):
        if letters[k] >= UInt8(ord("a")) and letters[k] <= UInt8(ord("z")):
            letters[k] -= 32
        if letters[k] == UInt8(ord(".")):
            letters[k] = GAP
    return letters^


def _rows(path: String, seqs: List[List[UInt8]]) raises -> List[List[UInt8]]:
    """Each distinct sequence's row of the MSA at `path`."""
    var msa = read_fasta(path)
    var rows = List[List[UInt8]]()
    for _ in range(len(seqs)):
        rows.append(List[UInt8]())
    for r in range(len(msa)):
        var row = _upper(List[UInt8](msa.seq(r)))
        var spelled = List[UInt8]()
        for x in row:
            if x != GAP:
                spelled.append(x)
        for d in range(len(seqs)):
            if len(rows[d]) == 0 and seqs[d] == spelled:
                rows[d] = row.copy()
    for d in range(len(seqs)):
        if len(rows[d]) == 0:
            raise Error(t"{path}: no row spells sequence {d}")
    return rows^


comptime INF = Int32(1 << 29)


def _partners(
    a: List[UInt8], b: List[UInt8], affine: Bool, mismatch: Int32, opening: Int32, extension: Int32
) -> Tuple[List[Int], List[Int]]:
    """For each letter of a, the least and most letter of b an optimal alignment of the pair matches it with (-1:
    none): at unit costs, or gap-affine (a gap of k letters `opening + k extension`, `mismatch` per substitution;
    Gotoh's three states forward, the match state's costs kept whole, and backward)."""
    var n = len(a)
    var m = len(b)
    var w = m + 1
    var lo_j = List[Int](length=n, fill=-1)
    var hi_j = List[Int](length=n, fill=-1)
    var mis = mismatch if affine else Int32(1)
    var oe = opening + extension if affine else Int32(1)
    var e = extension if affine else Int32(1)
    # Forward: FM whole; FX (a's letter against a gap) and FY (b's) in two rows. At unit costs the three states
    # collapse, the costs equal with opening 0 and extension 1.
    var fm = List[Int32](length=(n + 1) * w, fill=INF)
    var fx_prev = List[Int32](length=w, fill=INF)
    var fy_prev = List[Int32](length=w, fill=INF)
    var fx = List[Int32](length=w, fill=INF)
    var fy = List[Int32](length=w, fill=INF)
    fm[0] = 0
    for j in range(1, w):
        fy_prev[j] = oe + (Int32(j) - 1) * e
    for i in range(1, n + 1):
        fx[0] = oe + (Int32(i) - 1) * e
        fy[0] = INF
        fm[i * w] = INF
        for j in range(1, w):
            var sub = Int32(0) if a[i - 1] == b[j - 1] else mis
            var dm = fm[(i - 1) * w + j - 1]
            if fx_prev[j - 1] < dm:
                dm = fx_prev[j - 1]
            if fy_prev[j - 1] < dm:
                dm = fy_prev[j - 1]
            fm[i * w + j] = dm + sub if dm < INF else INF
            var gx = min(fm[(i - 1) * w + j] + oe, min(fx_prev[j] + e, fy_prev[j] + oe))
            fx[j] = gx
            var gy = min(fm[i * w + j - 1] + oe, min(fy[j - 1] + e, fx[j - 1] + oe))
            fy[j] = gy
        fx_prev = fx.copy()
        fy_prev = fy.copy()
    var best = min(fm[n * w + m], min(fx_prev[m], fy_prev[m]))
    # Backward, rolling: BM, BX, BY of row i from row i + 1.
    var bm_next = List[Int32](length=w + 1, fill=INF)
    var bx_next = List[Int32](length=w + 1, fill=INF)
    var by_next = List[Int32](length=w + 1, fill=INF)
    var bm = List[Int32](length=w + 1, fill=INF)
    var bx = List[Int32](length=w + 1, fill=INF)
    var by = List[Int32](length=w + 1, fill=INF)
    for i in reversed(range(n + 1)):
        for j in reversed(range(m + 1)):
            if i == n and j == m:
                bm[j] = 0
                bx[j] = 0
                by[j] = 0
                continue
            var to_m = INF
            if i < n and j < m:
                to_m = (Int32(0) if a[i] == b[j] else mis) + bm_next[j + 1]
            var to_x = bx_next[j] if i < n else INF
            var to_y = by[j + 1] if j < m else INF
            bm[j] = min(to_m, min(to_x + oe if to_x < INF else INF, to_y + oe if to_y < INF else INF))
            bx[j] = min(to_m, min(to_x + e if to_x < INF else INF, to_y + oe if to_y < INF else INF))
            by[j] = min(to_m, min(to_x + oe if to_x < INF else INF, to_y + e if to_y < INF else INF))
        # Letter a[i] (row i + 1) matched with b[j]: FM[i + 1][j + 1] + BM[i + 1][j + 1] == best.
        if i < n:
            for j in range(m):
                if fm[(i + 1) * w + j + 1] + bm_next[j + 1] == best:
                    if lo_j[i] < 0:
                        lo_j[i] = j
                    hi_j[i] = j
        bm_next = bm.copy()
        bx_next = bx.copy()
        by_next = by.copy()
    return (lo_j^, hi_j^)


def main() raises:
    var raw = argv()
    var args = List[String]()
    var affine = False
    for k in range(len(raw)):
        if String(raw[k]) == "--affine":
            affine = True
        else:
            args.append(String(raw[k]))
    if len(args) < 5:
        print("usage: structure <sequences.fa> <msa.fa> <period> <max pairs> [<msa.fa> ...]")
        return
    var records = read_fasta(String(args[1]))
    var period = Int(String(args[3]))
    var only_a = -1
    var only_b = -1
    var max_pairs = 0
    if String(args[4]) == "--pair":
        only_a = Int(String(args[5]))
        only_b = Int(String(args[6]))
        max_pairs = 1
    else:
        max_pairs = Int(String(args[4]))
    var seqs = List[List[UInt8]]()
    for r in range(len(records)):
        var letters = _upper(List[UInt8](records.seq(r)))
        var dup = False
        for x in seqs:
            if x == letters:
                dup = True
                break
        if not dup:
            seqs.append(letters^)
    var n_seqs = len(seqs)
    var paths = List[String]()
    paths.append(String(args[2]))
    for k in range(5 if only_a < 0 else 7, len(args)):
        paths.append(String(args[k]))
    var msas = List[List[List[UInt8]]]()
    for p in paths:
        msas.append(_rows(p, seqs))
    # Copies: the region's unit from the longest sequence, every array cut against it.
    var longest = 0
    for k in range(1, n_seqs):
        if len(seqs[k]) > len(seqs[longest]):
            longest = k
    var top = repeat_at(seqs[longest], period)
    if top.period == 0:
        print("no array at period", period)
        return
    var unit = consensus_unit(seqs[longest], top.start, top.end, period)
    var cut = List[Copies]()
    var arrays_lo = List[Int]()
    var arrays_hi = List[Int]()
    # For each sequence, each letter's copy (-1 outside the array).
    var copy_of = List[List[Int]]()
    var symbol = Dict[String, Int]()
    # Each distinct copy string's id, in order of first sighting (the copy map's).
    var ids = Dict[String, Int]()
    for k in range(n_seqs):
        var r = repeat_at(seqs[k], period)
        var cs = copies(seqs[k], r, unit)
        var which = List[Int](length=len(seqs[k]), fill=-1)
        for c in range(len(cs.starts)):
            for t in range(cs.starts[c], cs.ends[c]):
                which[t] = c
            var key = String(unsafe_from_utf8=List[UInt8](Span(seqs[k])[cs.starts[c] : cs.ends[c]]))
            symbol[key] = symbol.get(key, 0) + 1
            if key not in ids:
                ids[key] = len(ids)
        arrays_lo.append(r.start)
        arrays_hi.append(r.end)
        cut.append(cs^)
        copy_of.append(which^)
    var pairs = List[Tuple[Int, Int]]()
    for x in range(n_seqs):
        for y in range(x + 1, n_seqs):
            pairs.append((x, y) if len(seqs[x]) >= len(seqs[y]) else (y, x))
    var chosen = List[Tuple[Int, Int]]()
    if only_a >= 0:
        pairs.clear()
        pairs.append((only_a, only_b))
    var step = max(1, len(pairs) // max(1, max_pairs))
    var t = 0
    while t < len(pairs) and len(chosen) < max_pairs:
        chosen.append(pairs[t])
        t += step
    # Tallies per MSA: classes x (zone 0..2) x (rare 0/1); and right/total per sequence; flank letters.
    comptime CLASSES = 5  # right, off 1, off 2+, gapped, paired
    var n_msa = len(msas)
    var tally = List[Int](length=n_msa * CLASSES * 3 * 2, fill=0)
    var seq_right = List[Int](length=n_msa * n_seqs, fill=0)
    var seq_total = List[Int](length=n_msa * n_seqs, fill=0)
    var flank_on = List[Int](length=n_msa, fill=0)
    var flank_all = List[Int](length=n_msa, fill=0)
    # Shared unique copies: a copy string once in each of the pair's arrays (a copy-level MUM). How many, how many an
    # optimal alignment can pair, and per MSA how many it pairs.
    var shared_total = 0
    var shared_optimal = 0
    var shared_kept = List[Int](length=n_msa, fill=0)
    var done = 0
    for pq in chosen:
        var ia = pq[0]
        var ib = pq[1]
        ref a = seqs[ia]
        ref b = seqs[ib]
        var n = len(a)
        var m = len(b)
        if (n + 1) * (m + 1) > 150_000_000:
            continue
        done += 1
        var partners = _partners(a, b, affine, 4, 6, 2)
        var lo_j = partners[0].copy()
        var hi_j = partners[1].copy()
        ref ca = cut[ia]
        # The pair's copy-level MUMs: copy c of a and copy d of b, their string once in each.
        var count_a = Dict[String, Int]()
        var count_b = Dict[String, Int]()
        var where_b = Dict[String, Int]()
        for c in range(len(ca.starts)):
            var key = String(unsafe_from_utf8=List[UInt8](Span(a)[ca.starts[c] : ca.ends[c]]))
            count_a[key] = count_a.get(key, 0) + 1
        for d in range(len(cut[ib].starts)):
            var key = String(unsafe_from_utf8=List[UInt8](Span(b)[cut[ib].starts[d] : cut[ib].ends[d]]))
            count_b[key] = count_b.get(key, 0) + 1
            where_b[key] = d
        var mum_a = List[Int]()
        var mum_b = List[Int]()
        for c in range(len(ca.starts)):
            var key = String(unsafe_from_utf8=List[UInt8](Span(a)[ca.starts[c] : ca.ends[c]]))
            if count_a.get(key, 0) == 1 and count_b.get(key, 0) == 1:
                mum_a.append(c)
                mum_b.append(where_b[key])
        shared_total += len(mum_a)
        for u in range(len(mum_a)):
            var lo = Int.MAX
            var hi = -1
            for k in range(ca.starts[mum_a[u]], ca.ends[mum_a[u]]):
                if lo_j[k] >= 0:
                    for jj in [lo_j[k], hi_j[k]]:
                        var cbk = copy_of[ib][jj]
                        if cbk >= 0:
                            lo = min(lo, cbk)
                            hi = max(hi, cbk)
            if mum_b[u] >= lo and mum_b[u] <= hi:
                shared_optimal += 1
        for mi in range(n_msa):
            ref ra = msas[mi][ia]
            ref rb = msas[mi][ib]
            # Each letter of a: the letter of b in its column (-1: a gap there).
            var facing = List[Int](length=n, fill=-1)
            var i = 0
            var j = 0
            for c in range(len(ra)):
                var xa = ra[c] != GAP
                var xb = rb[c] != GAP
                if xa and xb:
                    facing[i] = j
                if xa:
                    i += 1
                if xb:
                    j += 1
            # Flanks: letters outside a's array, facing a letter of b an optimal alignment pairs them with.
            for k in range(n):
                if copy_of[ia][k] < 0:
                    flank_all[mi] += 1
                    var f = facing[k]
                    if f >= 0 and lo_j[k] >= 0 and f >= lo_j[k] and f <= hi_j[k]:
                        flank_on[mi] += 1
            # The shared unique copies the MSA pairs (most of a's copy's letters facing b's).
            for u in range(len(mum_a)):
                var c0 = ca.starts[mum_a[u]]
                var c1 = ca.ends[mum_a[u]]
                var hits = 0
                for k in range(c0, c1):
                    if facing[k] >= 0 and copy_of[ib][facing[k]] == mum_b[u]:
                        hits += 1
                if hits * 2 > c1 - c0:
                    shared_kept[mi] += 1
            for c in range(len(ca.starts)):
                var s0 = ca.starts[c]
                var s1 = ca.ends[c]
                # The optimal range of b's copies, and whether an optimal alignment pairs any of its letters.
                var range_lo = Int.MAX
                var range_hi = -1
                for k in range(s0, s1):
                    if lo_j[k] >= 0:
                        for jj in [lo_j[k], hi_j[k]]:
                            var cbk = copy_of[ib][jj]
                            if cbk >= 0:
                                range_lo = min(range_lo, cbk)
                                range_hi = max(range_hi, cbk)
                var pairable = range_hi >= 0
                # The MSA's partner: b's copy most of its letters face.
                var votes = Dict[Int, Int]()
                var facing_letters = 0
                for k in range(s0, s1):
                    var f = facing[k]
                    if f >= 0:
                        facing_letters += 1
                        var cbk = copy_of[ib][f]
                        if cbk >= 0:
                            votes[cbk] = votes.get(cbk, 0) + 1
                var partner = -1
                var partner_votes = 0
                for e in votes.items():
                    if e.value > partner_votes:
                        partner_votes = e.value
                        partner = e.key
                var cls = 0
                if facing_letters * 2 < s1 - s0 or partner < 0:
                    cls = 3 if pairable else 0
                elif not pairable:
                    cls = 4
                elif partner >= range_lo and partner <= range_hi:
                    cls = 0
                else:
                    var off = range_lo - partner if partner < range_lo else partner - range_hi
                    cls = 1 if off == 1 else 2
                var count = len(ca.starts)
                if only_a >= 0:
                    # The copy map: a's copy, its string's id, the optimal range in b, the MSA's partner and its id.
                    var key_a = String(unsafe_from_utf8=List[UInt8](Span(a)[s0:s1]))
                    var id_b = String("-")
                    if partner >= 0:
                        var pk = String(
                            unsafe_from_utf8=List[UInt8](Span(b)[cut[ib].starts[partner] : cut[ib].ends[partner]])
                        )
                        id_b = String(ids.get(pk, -1))
                    var verdict = ["right", "off1", "OFF", "GAPPED", "paired"][cls]
                    print(
                        "  copy",
                        c,
                        "id",
                        ids.get(key_a, -1),
                        "optimal",
                        String(range_lo) + "-" + String(range_hi) if pairable else String("gap"),
                        "msa",
                        partner,
                        "id",
                        id_b,
                        "facing",
                        facing_letters,
                        "of",
                        s1 - s0,
                        verdict,
                    )
                var zone = 0 if c * 10 < count else (2 if c * 10 >= 9 * count else 1)
                var key = String(unsafe_from_utf8=List[UInt8](Span(a)[s0:s1]))
                var rare = 1 if symbol.get(key, 0) <= 2 else 0
                tally[((mi * CLASSES + cls) * 3 + zone) * 2 + rare] += 1
                seq_total[mi * n_seqs + ia] += 1
                if cls == 0:
                    seq_right[mi * n_seqs + ia] += 1
    print("pairs", done, "of", len(pairs), "period", period, "unit", String(unsafe_from_utf8=unit.copy()))
    print(
        "shared unique copies (a copy string once in each of a pair):",
        shared_total,
        "| an optimal alignment can pair them:",
        shared_optimal,
    )
    var names = ["right", "off 1", "off 2+", "gapped", "paired"]
    for mi in range(n_msa):
        print("==", paths[mi])
        var total = 0
        for v in range(CLASSES * 6):
            total += tally[mi * CLASSES * 6 + v]
        for cls in range(CLASSES):
            var line = String(names[cls]) + ":"
            var all = 0
            for zone in range(3):
                for rare in range(2):
                    all += tally[((mi * CLASSES + cls) * 3 + zone) * 2 + rare]
            line += " " + String(all) + " (" + String(Float64(all * 1000 // max(1, total)) / 10.0) + "%)"
            for zone in range(3):
                var z = 0
                for rare in range(2):
                    z += tally[((mi * CLASSES + cls) * 3 + zone) * 2 + rare]
                line += " | " + ["start", "middle", "end"][zone] + " " + String(z)
            line += " | rare " + String(
                tally[((mi * CLASSES + cls) * 3 + 0) * 2 + 1]
                + tally[((mi * CLASSES + cls) * 3 + 1) * 2 + 1]
                + tally[((mi * CLASSES + cls) * 3 + 2) * 2 + 1]
            )
            print("  ", line)
        var rare_right = 0
        var rare_all = 0
        for cls in range(CLASSES):
            for zone in range(3):
                rare_all += tally[((mi * CLASSES + cls) * 3 + zone) * 2 + 1]
                if cls == 0:
                    rare_right += tally[((mi * CLASSES + cls) * 3 + zone) * 2 + 1]
        print(
            "   rare variant copies right:",
            rare_right,
            "of",
            rare_all,
            "| flank letters on an optimal alignment:",
            flank_on[mi],
            "of",
            flank_all[mi],
            "| shared unique copies paired:",
            shared_kept[mi],
            "of",
            shared_total,
        )
        var line = String("   copies right per sequence:")
        for k in range(n_seqs):
            if seq_total[mi * n_seqs + k] > 0:
                line += " " + String(k) + ":" + String(seq_right[mi * n_seqs + k] * 100 // seq_total[mi * n_seqs + k])
        print(line)
