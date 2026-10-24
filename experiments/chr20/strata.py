#!/usr/bin/env python3
"""Each structural error of an MSA's calls (truvari's fn and fp per region, from batch_call.py) put in a stratum:
one the MSA can fix, or one it cannot.

usage: python3 -I experiments/chr20/strata.py <batch_call out dir> <prep dir> <packages.jsonl.gz> <label>
           [--others a,b] [--regions file] [--jobs n] [--out strata.tsv]

Per region with an error, HG002's two truth haplotypes against the panel (every allele), the label's graph (its
sampled rows) and the label's two called haplotypes (haps.py); a haplotype within TOL edits (10 + 0.2% of its length)
counts as the same. Each truth haplotype, paired with the called one nearer it, is measured against the best the
graph allows (its row nearest the truth) and falls in the first stratum it meets:

  right       the call is as near the truth as the graph allows, and the truth is in the panel: the error is how
              the records split a copy-number change (one insertion of two copies against the truth's two of one),
              which truvari matches record by record: the haplotype is right
  best        the same, the truth in no panel allele: nothing an MSA of these alleles can call
  sampled     the call is worse, and the panel holds an allele nearer the truth than any the sampling kept: the
              graph's input, not its alignment
  mosaic      the call is worse and no row of the graph: the genotyper joined rows across snarls the MSA cut
  confident   the call is another row, every indel call of the region at GQ > 10: the graph misled the genotyper
  coinflip    the call is another row at GQ <= 10: reads that cannot tell the alleles apart

`~absent` marks a stratum past `best` whose truth is in no panel allele; with --others, `+other` marks one where
another label's call of the same reads is the truth (a graph can get it). A mosaic, confident or coinflip call is
also marked `reads` when the pair of the graph's rows the region's reads support best by K-mer counts (oracle.py:
depth as well as presence) is the truth, as near as the graph allows, and `noinfo` when it is not: the reads cannot
tell the alleles apart (often copies added to an array, which reads of the reference's copies fit as well).
Fixable on the MSA's side: `reads` and `+other`; not: right (representation), best, sampled, `noinfo`.
Each fn record is charged to the truth haplotypes carrying it (half to each when both do), each fp to the truth
haplotype its called haplotype pairs with. Prints errors per stratum and the regions, worst first, of the fixable
ones; writes one line per error to --out.
"""
import gzip
import hashlib
import json
import multiprocessing
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import haps  # noqa: E402
import oracle  # noqa: E402

FIXABLE = ("reads",)


def tol(seq):
    return 10 + len(seq) // 500


def gt_of(sample):
    return sample.split(":")[0].replace("/", "|").split("|")


def errors(out_dir, label, rid):
    """(kind, pos, ref length, alt length, genotype) of truvari's fn and fp records of a region's structural calls."""
    out = []
    for kind in ("fn", "fp"):
        path = os.path.join(out_dir, "%s.%s" % (label, rid), "score", "chr20.truvari", kind + ".vcf.gz")
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                w = line.rstrip("\n").split("\t")
                out.append((kind, int(w[1]), len(w[3]), max(len(a) for a in w[4].split(",")), gt_of(w[9])))
    return out


def min_gq(calls_path):
    gq = []
    for line in open(calls_path):
        if line.startswith("#"):
            continue
        w = line.rstrip("\n").split("\t")
        if abs(max(len(a) for a in w[4].split(",")) - len(w[3])) < 10:
            continue
        f = dict(zip(w[8].split(":"), w[9].split(":")))
        gq.append(int(f["GQ"]) if f.get("GQ", ".").lstrip("-").isdigit() else 0)
    return min(gq) if gq else 99


def stratify(task):
    out_dir, prep_dir, pkg, label, others = task
    rid = pkg["id"]
    errs = errors(out_dir, label, rid)
    if not errs:
        return rid, [], None
    try:
        meta, contig, truth = haps.region(out_dir, prep_dir, pkg, rid)
        weight = {a["id"]: a["weight"] for a in pkg["alleles"]}
        panel = {}
        for a in pkg["alleles"]:
            panel.setdefault(meta["left"] + a["seq"].upper() + meta["right"], []).append(a["id"])
        calls, sampled, called, pair = haps.called(out_dir, rid, label, meta, contig, truth)
        rows = json.load(open(os.path.join(out_dir, "%s.%s" % (label, rid), "rows.json")))
        key = hashlib.sha1(json.dumps([rid, sorted(rows.items())]).encode()).hexdigest()[:20]
        gq = min_gq(os.path.join(out_dir, "cache", key + ".vcf"))
        other_right = [False, False]
        for o in others:
            try:
                _, _, ocalled, opair = haps.called(out_dir, rid, o, meta, contig, truth)
            except (OSError, ValueError, KeyError):
                continue
            for h in (0, 1):
                t = opair[h]
                if haps.distance(ocalled[h], truth[t]) <= tol(truth[t]):
                    other_right[t] = True
        # The pair of the graph's rows the reads' K-mer counts support best, paired with the truth.
        oa, ob, margin, sa, sb = oracle.best_pair(rows, meta["left"], meta["right"], prep_dir, rid)
        straight = haps.distance(sa, truth[0]) + haps.distance(sb, truth[1])
        crossed = haps.distance(sa, truth[1]) + haps.distance(sb, truth[0])
        reads_pick = [sa, sb] if straight <= crossed else [sb, sa]
        stratum = {}
        info = {}
        for h in (0, 1):
            t = pair[h]
            T, C = truth[t], called[h]
            d_called = haps.distance(C, T)
            d_panel, a_panel, _ = haps.nearest(T, panel, weight)
            d_graph, a_graph, _ = haps.nearest(T, sampled, weight)
            d_row, a_row, _ = haps.nearest(C, sampled, weight)
            # Against the best the graph allows (its row nearest the truth): at it, the call is all the MSA's
            # alleles can give; past it, why.
            if d_called <= d_graph + tol(T):
                s = "right" if d_panel <= tol(T) else "best"
            elif d_graph > d_panel + tol(T):
                s = "sampled"
            elif d_row > tol(C):
                s = "mosaic"
            elif gq > 10:
                s = "confident"
            else:
                s = "coinflip"
            if s.split()[0] in ("mosaic", "confident", "coinflip"):
                # Do the reads hold the answer? Their best pair of rows is as near the truth as the graph allows.
                s += " reads" if haps.distance(reads_pick[t], T) <= d_graph + tol(T) else " noinfo"
            if d_panel > tol(T) and s not in ("best",):
                s += " ~absent"
            if other_right[t] and s.split()[0] not in ("right", "best"):
                s += " +other"
            stratum[t] = s
            info[t] = (d_called, d_panel, a_panel, d_graph, a_graph, d_row, a_row)
        out = []
        for kind, pos, rl, al, gt in errs:
            if kind == "fn":
                carriers = [t for t in (0, 1) if t < len(gt) and gt[t] not in ("0", ".")] or [0, 1]
            else:
                carriers = [pair[h] for h in (0, 1) if h < len(gt) and gt[h] not in ("0", ".")] or [0, 1]
            for t in carriers:
                out.append((rid, label, kind, pos, al - rl, t, stratum[t], 1.0 / len(carriers)))
        return rid, out, {"gq": gq, "strata": stratum, "info": info}
    except (OSError, ValueError, KeyError) as error:
        return rid, [(rid, label, kind, pos, al - rl, -1, "error", 1.0) for kind, pos, rl, al, _ in errs], str(error)


def fixable(s):
    return any(x in s.split() for x in FIXABLE) or "+other" in s


def main():
    args = sys.argv[1:]

    def option(name, default):
        if name in args:
            i = args.index(name)
            value = args[i + 1]
            del args[i : i + 2]
            return value
        return default

    others = [x for x in option("--others", "").split(",") if x]
    regions_file = option("--regions", "")
    jobs = int(option("--jobs", "8"))
    out_path = option("--out", "")
    out_dir, prep_dir, packages, label = args[:4]
    regions = set(open(regions_file).read().split()) if regions_file else None
    pkgs = []
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if (regions is None or p["id"] in regions) and os.path.isdir(os.path.join(out_dir, "%s.%s" % (label, p["id"]))):
                pkgs.append(p)
    tasks = [(out_dir, prep_dir, p, label, others) for p in sorted(pkgs, key=lambda p: -len(p["alleles"]))]
    rows, notes = [], {}
    with multiprocessing.Pool(jobs) as pool:
        for rid, out, note in pool.imap_unordered(stratify, tasks):
            rows += out
            if note is not None:
                notes[rid] = note
    totals = {}
    for r in rows:
        x = totals.setdefault(r[6], [0.0, 0.0])
        x[0 if r[2] == "fn" else 1] += r[7]
    print("%s: %d structural errors in %d regions (fn %.0f, fp %.0f)" % (
        label, len(set((r[0], r[2], r[3]) for r in rows)), len(set(r[0] for r in rows)),
        sum(x[0] for x in totals.values()), sum(x[1] for x in totals.values())))
    print("  %-22s %7s %7s %7s" % ("stratum", "fn", "fp", "all"))
    for s in sorted(totals, key=lambda s: (not fixable(s), -sum(totals[s]))):
        print("  %-22s %7.1f %7.1f %7.1f%s" % (s, totals[s][0], totals[s][1], sum(totals[s]),
                                              "  fixable" if fixable(s) else ""))
    fix = sum(sum(v) for s, v in totals.items() if fixable(s))
    print("  fixable on the MSA's side %.1f of %.1f" % (fix, sum(sum(v) for v in totals.values())))
    by_region = {}
    for r in rows:
        if fixable(r[6]):
            by_region[r[0]] = by_region.get(r[0], 0) + r[7]
    print("fixable regions, worst first:")
    for rid in sorted(by_region, key=lambda x: -by_region[x])[:25]:
        n = notes.get(rid)
        desc = ""
        if isinstance(n, dict):
            desc = "GQmin %d | " % n["gq"] + " | ".join(
                "t%d %s: called %d, panel %d (%s), graph %d (%s), call's row %d (%s)" % ((t, n["strata"][t]) + n["info"][t])
                for t in (0, 1))
        print("  %-20s %5.1f  %s" % (rid, by_region[rid], desc))
    if out_path:
        with open(out_path, "w") as f:
            f.write("region\tlabel\tkind\tpos\tsvlen\ttruth_hap\tstratum\tshare\n")
            for r in sorted(rows):
                f.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
