#!/usr/bin/env python3
"""Where an MSA's implied pairwise alignments lose to the optimal ones, per region, pair and stretch.

usage: python3 -I experiments/vntreval/pairs.py <vntr-eval> <karan> <msa dir> <regions> [--jobs n] [--out tsv]

<msa dir> holds `<id>.msa.fa` (rows named as hap32.fa's); <regions> is a comma list of region ids or a stratum
name. For every pair of distinct haplotypes in a region: the unit cost of the alignment the MSA implies (each
column where the two rows differ costs one, a column where both are gaps nothing), the optimal unit cost (karan),
and the implied cost's excess split by where it falls in CHM13's row: before the repeat's core, in it, or after
it. Prints per region the summed implied cost over the optimum (as vntr-eval's all_cost_over_opt counts it,
distinct pairs once), the share of the excess in the core, and how the excess spreads over pairs (the share of
pairs at the optimum, the worst pair). `--out` writes every pair's row.
"""
import csv
import itertools
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool


def read_rows(path):
    names, rows = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            rows.append([])
        elif line:
            rows[-1].append(line.upper().replace(".", "-"))
    return names, ["".join(r) for r in rows]


def region_profile(task):
    vntr_eval, karan, msa_dir, rid = task
    meta = json.load(open(os.path.join(vntr_eval, "regions", rid, "region.json")))
    path = os.path.join(msa_dir, rid + ".msa.fa")
    if not os.path.exists(path):
        return rid, None
    names, rows = read_rows(path)
    # The repeat core in MSA columns, from CHM13's row.
    chm13 = next((k for k, n in enumerate(names) if n.upper().startswith("CHM13")), None)
    core = None
    if chm13 is not None:
        offset_start = meta.get("core_start", 0) - meta.get("span_start", 0)
        offset_end = meta.get("core_end", 0) - meta.get("span_start", 0)
        letters = 0
        lo = hi = None
        for c, x in enumerate(rows[chm13]):
            if x != "-":
                if letters == offset_start and lo is None:
                    lo = c
                if letters == offset_end:
                    hi = c
                letters += 1
        core = (lo or 0, hi if hi is not None else len(rows[chm13]))
    seqs = [r.replace("-", "") for r in rows]
    distinct = {}
    for k, s in enumerate(seqs):
        distinct.setdefault(s, k)
    keep = sorted(distinct.values())
    pairs = list(itertools.combinations(keep, 2))
    with tempfile.TemporaryDirectory() as tmp:
        tsv = os.path.join(tmp, "pairs.tsv")
        with open(tsv, "w") as f:
            for n, (a, b) in enumerate(pairs):
                f.write("%d\t%s\t%s\n" % (n, seqs[a] or "A", seqs[b] or "A"))
        out = subprocess.run([karan, "--pairs", tsv, "--distance"], capture_output=True, text=True).stdout
    optimum = {}
    for line in out.strip().splitlines()[1:]:
        f = line.split("\t")
        optimum[int(f[0])] = int(f[2])
    result = []
    for n, (a, b) in enumerate(pairs):
        ra, rb = rows[a], rows[b]
        cost = [0, 0, 0]
        for c, (x, y) in enumerate(zip(ra, rb)):
            if x != y:
                where = 1 if core and core[0] <= c < core[1] else (0 if core and c < core[0] else 2)
                cost[where] += 1
        opt = optimum.get(n, max(len(seqs[a]), len(seqs[b])) if not seqs[a] or not seqs[b] else 0)
        result.append((names[a], names[b], len(seqs[a]), len(seqs[b]), opt, sum(cost), cost[0], cost[1], cost[2]))
    return rid, (meta.get("stratum"), meta.get("period"), result)


def main():
    args = sys.argv[1:]
    jobs = int(args.pop(args.index("--jobs") + 1)) if "--jobs" in args else 8
    if "--jobs" in args:
        args.remove("--jobs")
    out_path = args.pop(args.index("--out") + 1) if "--out" in args else None
    if "--out" in args:
        args.remove("--out")
    vntr_eval, karan, msa_dir, which = args
    if "," in which or which.startswith("L") or which.startswith("TR"):
        regions = which.split(",")
    else:
        with open(os.path.join(vntr_eval, "regions", "regions.tsv")) as f:
            regions = [r["region_id"] for r in csv.DictReader(f, delimiter="\t") if r["stratum"] == which]
    with Pool(jobs) as pool:
        profiles = pool.map(region_profile, [(vntr_eval, karan, msa_dir, rid) for rid in regions])
    writer = open(out_path, "w") if out_path else None
    if writer:
        writer.write("region\tstratum\ta\tb\ta_len\tb_len\toptimum\timplied\tleft\tcore\tright\n")
    print("| region | stratum | period | pairs | implied / optimum | excess in core | pairs at optimum | worst pair |")
    print("| :-- | :-- | --: | --: | --: | --: | --: | --: |")
    total_opt = total_imp = total_core = total_excess = 0
    for rid, profile in profiles:
        if profile is None:
            print("| %s | missing |" % rid)
            continue
        stratum, period, rows = profile
        opt = sum(r[4] for r in rows)
        imp = sum(r[5] for r in rows)
        excess = imp - opt
        core = sum(r[7] for r in rows) - 0
        # The core's share of the excess: the core's implied cost less its share of the optimum, approximated by
        # the core's share of the implied cost when the pair is at its optimum.
        core_excess = sum(max(0, r[7] - (r[4] * r[7] / r[5] if r[5] else 0)) for r in rows if r[5] > r[4])
        at_opt = sum(r[5] <= r[4] for r in rows)
        worst = max((r[5] / r[4] if r[4] else 1.0) for r in rows) if rows else 1.0
        total_opt += opt
        total_imp += imp
        total_core += core_excess
        total_excess += excess
        print("| %s | %s | %s | %d | %.4f | %.2f | %.2f | %.3f |" % (
            rid, stratum, period, len(rows), imp / opt if opt else 1.0,
            core_excess / excess if excess > 0 else 0.0, at_opt / len(rows) if rows else 1.0, worst))
        if writer:
            for r in rows:
                writer.write("%s\t%s\t%s\n" % (rid, stratum, "\t".join(str(x) for x in r)))
    print("\nall: implied / optimum %.4f; excess in the core %.2f" % (
        total_imp / total_opt if total_opt else 1.0, total_core / total_excess if total_excess > 0 else 0.0))


if __name__ == "__main__":
    main()
