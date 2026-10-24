#!/usr/bin/env python3
"""Which tandem-repeat regions an arm's structural calls get wrong, and how, against the truth and another arm
(truvari's per-variant output from experiments/chr20/score_calls.py).

usage: python3 -I experiments/chr20/call_sites.py <work dir> <packages.jsonl.gz> <arm> <other arm> [<set>] [--top n]

For every truth variant (truvari's tp-base, fn) and call (tp-comp, fp) of 50 bp or more in <set> (`repeats`), its
region (the package whose CHM13 span holds it), its type (INS or DEL by length change) and size. Per region and arm:
truth found (TP), truth missed (FN), calls not in the truth (FP), with the calls' sizes. Prints the totals, the
regions where <arm> errs and <other arm> does not, and the reverse, worst first, each with its truth variants and
both arms' calls (type, size, position) so a region can be looked at; writes <work dir>/call_sites.<arm>.<other>.tsv.
"""
import bisect
import collections
import gzip
import json
import os
import sys


def regions_of(packages):
    out = []
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p.get("contig") == "chr20" and p.get("ref_end"):
                out.append((int(p["ref_start"]), int(p["ref_end"]), p["id"], p.get("classes", "")))
    out.sort()
    return out


def records(path):
    """(position, type, size, length change) of each record of a VCF (gzipped), the longest alternative's."""
    out = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            w = line.split("\t", 6)
            pos, ref, alts = int(w[1]), w[3], w[4].split(",")
            best = max(alts, key=lambda a: abs(len(a) - len(ref)))
            change = len(best) - len(ref)
            kind = "INS" if change > 0 else "DEL" if change < 0 else "SUB"
            out.append((pos, kind, max(abs(change), len(ref) if kind == "SUB" else 0), change))
    return out


def place(regions, pos):
    starts = [r[0] for r in regions]
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0 and pos <= regions[i][1]:
        return regions[i][2]
    return None


def main():
    args = sys.argv[1:]
    top = 15
    if "--top" in args:
        top = int(args[args.index("--top") + 1])
        del args[args.index("--top") : args.index("--top") + 2]
    work, packages, arm, other = args[:4]
    which = args[4] if len(args) > 4 else "repeats"
    regions = regions_of(packages)
    starts = [r[0] for r in regions]
    classes = {r[2]: r[3] for r in regions}

    def where(pos):
        i = bisect.bisect_right(starts, pos) - 1
        return regions[i][2] if i >= 0 and pos <= regions[i][1] else None

    found = {}
    for a in (arm, other):
        tdir = os.path.join(work, "calls", "%s.%s" % (a, which), "score", "chr20.truvari")
        per = collections.defaultdict(lambda: collections.defaultdict(list))
        for part in ("tp-base", "fn", "tp-comp", "fp"):
            for rec in records(os.path.join(tdir, part + ".vcf.gz")):
                rid = where(rec[0])
                if rid is not None:
                    per[rid][part].append(rec)
        found[a] = per
    totals = {a: collections.Counter() for a in (arm, other)}
    for a in (arm, other):
        for rid, parts in found[a].items():
            for part, recs in parts.items():
                totals[a][part] += len(recs)
                for r in recs:
                    totals[a][part + ":" + r[1]] += 1
    for a in (arm, other):
        t = totals[a]
        print("%-9s truth found %d (INS %d, DEL %d), missed %d (INS %d, DEL %d); calls not in the truth %d (INS %d, "
              "DEL %d)" % (a, t["tp-base"], t["tp-base:INS"], t["tp-base:DEL"], t["fn"], t["fn:INS"], t["fn:DEL"],
                           t["fp"], t["fp:INS"], t["fp:DEL"]))
    rows = []
    for rid in set(found[arm]) | set(found[other]):
        pa, po = found[arm].get(rid, {}), found[other].get(rid, {})
        err_a = len(pa.get("fn", [])) + len(pa.get("fp", []))
        err_o = len(po.get("fn", [])) + len(po.get("fp", []))
        rows.append((err_a - err_o, rid, pa, po, err_a, err_o))
    with open(os.path.join(work, "call_sites.%s.%s.tsv" % (arm, other)), "w") as f:
        f.write("region\tclasses\t%s_errors\t%s_errors\t%s_fn\t%s_fp\t%s_fn\t%s_fp\n" % (arm, other, arm, arm, other,
                                                                                         other))
        for d, rid, pa, po, ea, eo in sorted(rows, reverse=True):
            f.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (rid, classes.get(rid, ""), ea, eo, len(pa.get("fn", [])),
                                                          len(pa.get("fp", [])), len(po.get("fn", [])),
                                                          len(po.get("fp", []))))
    worse = sum(1 for r in rows if r[0] > 0)
    better = sum(1 for r in rows if r[0] < 0)
    print("regions with a structural error in either: %d; %s errs more in %d, %s in %d" % (len(rows), arm, worse,
                                                                                         other, better))

    def show(recs):
        return " ".join("%s%d@%d" % (r[1], r[2], r[0]) for r in recs[:6]) + (" ..." if len(recs) > 6 else "")

    rows.sort(reverse=True)
    for title, chosen in (("where %s errs and %s does less" % (arm, other), rows[:top]),
                          ("where %s errs and %s does less" % (other, arm), rows[::-1][:top])):
        print(title + ":")
        for d, rid, pa, po, ea, eo in chosen:
            if d == 0:
                continue
            truth = sorted(pa.get("tp-base", []) + pa.get("fn", []))
            print("  %-20s %-10s %s errors %d (missed %d, spurious %d), %s %d | truth: %s" % (
                rid, classes.get(rid, "")[:10], arm, ea, len(pa.get("fn", [])), len(pa.get("fp", [])), other, eo,
                show(truth)))
            print("  %20s %s calls: %s | %s calls: %s" % (
                "", arm, show(sorted(pa.get("tp-comp", []) + pa.get("fp", []))), other,
                show(sorted(po.get("tp-comp", []) + po.get("fp", [])))))


if __name__ == "__main__":
    main()
