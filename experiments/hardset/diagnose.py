"""Which sites an MSA gets wrong on a hard-set region, and why (experiments/hardset/hardset.py --diagnose).

For each pair qc scores, every letter of its first allele is checked against qc's optimum (the reference alignment
pair F1 counts) and against every optimal alignment of the pair (build/anchors --coopt: two-piece forward and
backward DP under pgrealign's scores):
- `right`: aligned as the reference aligns it;
- `tie`: aligned otherwise but as some optimal alignment aligns it (an equally good placement the reference puts
  elsewhere: a tie-breaking convention would recover it);
- `wrong`: aligned as no optimal alignment aligns it, by kind against the reference partner: `near` (under 3 letters
  away), `copy` (a whole number of periods), `far`, or `gapped` (the reference leaves it unaligned).
Every letter is placed on the region's reference coordinate (CHM13's allele through the letter's MSA column, plus the
region's CHM13 start). A pair's consecutive wrong letters with one offset form a block, one misplaced event: a copy's
indel placed elsewhere shifts every letter between the optimum's place and the MSA's by a copy, and the block's ends
are those two places. Blocks of the same kind and shift that overlap, over all pairs, are a site: its CHM13 interval
(its id), wrong letters, ties within it, the pairs and alleles involved, its kind and shift, the star's wrong letters
on the same interval, and its longest block spelled out (where the MSA and the optimum send those letters).

Per region the summary: f1 (qc), co_f1 (any optimal alignment counted right), the wrong letters by kind, `top3` (their
share in the three alleles most involved), `excess` (the induced alignments' summed cost over the optima) and `sp`
(the MSA's summed cost under pgrealign's scores over up to SP_PAIRS weighted allele pairs, over the star's: under 1
with a lower f1 means the objective prefers this MSA and F1 does not; over 1 the search found a worse MSA by its own
measure too).
"""
import collections
import os
import random
import statistics
import subprocess
import tempfile

SP_PAIRS = 300


def partners(ra, rb):
    """Each letter of the first sequence: its partner in the second, or None."""
    out = []
    j = 0
    for x, y in zip(ra, rb):
        if x != "-":
            out.append(j if y != "-" else None)
        if y != "-":
            j += 1
    return out


def coopt(anchors, lines):
    """build/anchors --coopt on `lines` [(a, b, row_a, row_b)]: [(best, cost, pairs, on, flags)]."""
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        for a, b, ra, rb in lines:
            f.write("%s\t%s\t%s\t%s\n" % (a, b, ra, rb))
        path = f.name
    try:
        out = subprocess.run([anchors, "--coopt", path], capture_output=True, text=True, check=True).stdout
    finally:
        os.unlink(path)
    res = []
    for line in out.splitlines():
        w = line.split()
        res.append((int(w[0]), int(w[1]), int(w[2]), int(w[3]), w[4] if len(w) > 4 else ""))
    return res


def sp_cost(hs, qc, rows, pairs):
    return sum(hs.pair_cost(*qc.induced(rows[a], rows[b])) for a, b in pairs)


def sp_pairs(alleles, seed):
    """Up to SP_PAIRS pairs of distinct alleles, drawn with the alleles' weights (pairs of haplotypes)."""
    ids = [a for a, _, _ in alleles]
    weights = [w for _, _, w in alleles]
    if len(ids) < 2:
        return []
    rng = random.Random(seed)
    out = set()
    tries = 0
    while len(out) < min(SP_PAIRS, len(ids) * (len(ids) - 1) // 2) and tries < 20 * SP_PAIRS:
        a, b = rng.choices(ids, weights=weights, k=2)
        tries += 1
        if a != b:
            out.add((min(a, b), max(a, b)))
    return sorted(out)


def ref_frame(rows, ref, ref_start):
    """Per MSA column, the CHM13 coordinate there (the reference allele's next letter), and per allele its letters'
    columns."""
    width = len(next(iter(rows.values())))
    coord = []
    k = 0
    row = rows.get(ref)
    for c in range(width):
        coord.append(ref_start + k)
        if row is not None and row[c] != "-":
            k += 1
    columns = {a: [c for c, x in enumerate(r) if x != "-"] for a, r in rows.items()}
    return coord, columns


def wrong_letters(hs, qc, anchors, rows, opt, period, ref, ref_start):
    """Every scored pair's letters, judged: [(pair, coord, letter, state, kind, shift)] and the pair records."""
    lines, recs = [], []
    for rec in opt:
        if rec.get("skipped") or rec["kind"] != "pair":
            continue
        oa, ob = rec["rows"]
        ra, rb = qc.induced(rows[rec["a"]], rows[rec["b"]])
        lines.append((oa.replace("-", ""), ob.replace("-", ""), ra, rb))
        recs.append((rec, oa, ob, ra, rb))
    found = coopt(anchors, lines)
    coord, columns = ref_frame(rows, ref, ref_start)
    letters = []
    for (rec, oa, ob, ra, rb), (best, cost, pairs, on, flags) in zip(recs, found):
        ref_p, mine = partners(oa, ob), partners(ra, rb)
        cols = columns[rec["a"]]
        for i, flag in enumerate(flags):
            r, x = ref_p[i], mine[i]
            if flag == "." or r == x:
                continue
            state = "tie" if flag == "1" else "wrong"
            kind, shift = "", 0
            if r is None:
                kind = "gapped"
            else:
                off = x - r
                shift = round(off / period) if period else 0
                if abs(off) < 3:
                    kind = "near"
                elif period > 0 and abs(off) >= period - 2 and min(abs(off) % period, period - abs(off) % period) <= 2:
                    kind = "copy"
                else:
                    kind = "far"
            letters.append(((rec["a"], rec["b"]), coord[cols[i]], i, state, kind, shift, x, r))
    return letters, recs, found


def blocks_of(letters):
    """Per pair, the runs of consecutive wrong letters with one offset (the MSA's partner minus the optimum's; None
    for a letter the optimum leaves unaligned): one misplaced event each, e.g. a copy's indel placed elsewhere, the
    letters between the two placements shifted. [(pair, start coord, end coord, first letter, last letter, kind,
    shift, offset, MSA partners, optimum partners)]."""
    by_pair = collections.defaultdict(list)
    for x in letters:
        if x[3] == "wrong":
            by_pair[x[0]].append(x)
    out = []
    for pair, xs in by_pair.items():
        xs.sort(key=lambda x: x[2])
        run = [xs[0]]
        for x in xs[1:]:
            prev = run[-1]
            off = None if x[7] is None else x[6] - x[7]
            prev_off = None if prev[7] is None else prev[6] - prev[7]
            if x[2] == prev[2] + 1 and off == prev_off:
                run.append(x)
            else:
                out.append(_block(pair, run))
                run = [x]
        out.append(_block(pair, run))
    return out


def _block(pair, run):
    off = None if run[0][7] is None else run[0][6] - run[0][7]
    coords = [x[1] for x in run]
    return (pair, min(coords), max(coords), run[0][2], run[-1][2], run[0][4], run[0][5], off,
            (run[0][6], run[-1][6]), (run[0][7], run[-1][7]))


def sites_of(blocks):
    """Blocks of every pair merged into sites: same kind and shift, overlapping CHM13 intervals."""
    out = []
    for b in sorted(blocks, key=lambda b: (b[5], b[6], b[1])):
        if out and out[-1][-1][5] == b[5] and out[-1][-1][6] == b[6] and b[1] <= max(x[2] for x in out[-1]):
            out[-1].append(b)
        else:
            out.append([b])
    return out


def describe(site, ties, star_wrong, period):
    """A site's row: its interval, the letters and pairs in it, its kind and shift, the star's wrong letters there and
    its longest block spelled out."""
    start, end = min(b[1] for b in site), max(b[2] for b in site)
    letters = sum(b[4] - b[3] + 1 for b in site)
    pairs = collections.Counter(b[0] for b in site)
    alleles = collections.Counter()
    for (a, b), c in pairs.items():
        alleles[a] += c
        alleles[b] += c
    longest = max(site, key=lambda b: b[4] - b[3])
    (a, b), kind, shift = longest[0], longest[5], longest[6]
    if kind == "gapped":
        what = "the MSA aligns %s letters %d-%d to %s %d-%d, the optimum leaves them unaligned" % (
            a, longest[3], longest[4], b, longest[8][0], longest[8][1])
    else:
        what = "the MSA aligns %s letters %d-%d to %s %d-%d, the optimum to %d-%d (%+d letters%s)" % (
            a, longest[3], longest[4], b, longest[8][0], longest[8][1], longest[9][0], longest[9][1], longest[7],
            ", %+d copies" % shift if kind == "copy" else "")
    tie_here = sum(1 for x in ties if start <= x[1] <= end)
    star_here = sum(1 for p in star_wrong if start <= p <= end) if star_wrong is not None else None
    return {
        "start": start,
        "end": end,
        "wrong": letters,
        "ties": tie_here,
        "pairs": len(pairs),
        "alleles": ",".join("%s:%d" % kv for kv in alleles.most_common(3)),
        "kinds": kind,
        "shift": "%+d" % shift if kind == "copy" else "",
        "star_wrong": star_here,
        "worst": what,
        "pair_list": ";".join("%s~%s" % ab for ab in sorted(pairs)),
    }


def diagnose(hs, qc, anchors, rows, opt, period, alleles, ref, ref_start, star_rows=None):
    """The region's summary and sites for one MSA (`rows`: {allele id: row}); `star_rows`, the star's, for `sp` and
    the star's wrong letters on each site. Also returns the wrong letters' coordinates (for --against)."""
    letters, recs, found = wrong_letters(hs, qc, anchors, rows, opt, period, ref, ref_start)
    tp = n_o = n_m = on = excess = 0
    for (rec, oa, ob, ra, rb), (best, cost, pairs, on_opt, flags) in zip(recs, found):
        o, m = qc.aligned_pairs(oa, ob), qc.aligned_pairs(ra, rb)
        tp += len(o & m)
        n_o += len(o)
        n_m += len(m)
        on += on_opt
        excess += cost - best
    wrong = [x for x in letters if x[3] == "wrong"]
    ties = [x for x in letters if x[3] == "tie"]
    by_allele = collections.Counter()
    for x in wrong:
        by_allele[x[0][0]] += 1
        by_allele[x[0][1]] += 1
    kinds = collections.Counter(x[4] for x in wrong)
    summary = {
        "f1": 2 * tp / (n_o + n_m) if n_o + n_m else None,
        "co_f1": 2 * on / (n_o + n_m) if n_o + n_m else None,
        "wrong": len(wrong),
        "ties": len(ties),
        "near": kinds["near"],
        "copy": kinds["copy"],
        "far": kinds["far"],
        "gapped": kinds["gapped"],
        "top3": sum(c for _, c in by_allele.most_common(3)) / (2 * len(wrong)) if wrong else 0.0,
        "excess": excess,
    }
    star_wrong = None
    if star_rows is not None:
        pairs = sp_pairs(alleles, len(alleles))
        mine_sp, star_sp = sp_cost(hs, qc, rows, pairs), sp_cost(hs, qc, star_rows, pairs)
        summary["sp"] = mine_sp / star_sp if star_sp else None
        if star_rows is not rows:
            star_letters, _, _ = wrong_letters(hs, qc, anchors, star_rows, opt, period, ref, ref_start)
            star_wrong = [x[1] for x in star_letters if x[3] == "wrong"]
    sites = [describe(s, ties, star_wrong, period) for s in sites_of(blocks_of(letters))]
    sites.sort(key=lambda s: -s["wrong"])
    return summary, sites, sorted([x[0][0], x[0][1], x[1]] for x in wrong)


def against(sites_before, wrong_before, wrong_now):
    """Each earlier site re-scored: its pairs' wrong letters within its interval, before and now (`wrong_*`: [a, b,
    coord] per wrong letter); fixed (under a tenth left), better (under 0.8), same, or worse (over 1.25)."""

    def count(wrong, site, pairs):
        return sum(1 for a, b, c in wrong if (a, b) in pairs and site["start"] <= c <= site["end"])

    out = []
    for site in sites_before:
        pairs = {tuple(x.split("~")) for x in site["pair_list"].split(";") if x}
        before, now = count(wrong_before, site, pairs), count(wrong_now, site, pairs)
        status = ("fixed" if now <= before / 10 else "better" if now < 0.8 * before else
                  "worse" if now > 1.25 * before else "same")
        out.append((site, before, now, status))
    return out
