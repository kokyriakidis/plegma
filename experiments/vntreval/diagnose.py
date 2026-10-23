#!/usr/bin/env python3
"""Where a plegma method loses to another MSA on vntr-eval, and why: its search, or its objective.

usage: python3 -I experiments/vntreval/diagnose.py <vntr-eval> <method> <other msa dir> <other method> [--jobs n]

For every region with both MSAs (<vntr-eval>/work/plegma/<method>/<id>.msa.fa and <other msa dir>/<id>.msa.fa)
and both scores (results/<method>/, results/<other method>/): each MSA's sum of pairs over every pair of its rows (a
pair of distinct sequences weighing the product of their rows, as vntr-eval counts) at unit costs, plegma's objective (unit plus 1 per gap a pair opens) and vntr-eval's affine costs (sp_score.py), and
the two methods' all_cost_over_opt and affine_all_cost_over_opt. A region the other method wins on the unit metric
is
- `search` when the other MSA is lower on plegma's objective too: plegma's search missed a better MSA;
- `objective` when plegma's MSA is lower on its objective: the objective prefers what the metric does not;
and likewise for the affine metric. Prints the classes per stratum, and the regions lost by most.
"""
import collections
import csv
import json
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp_score import distinct, pair, read_rows  # noqa: E402


def costs(path):
    """The MSA's sums of pairs over every pair of its rows, as vntr-eval counts pairs: a pair of distinct sequences
    weighs the product of how many rows spell each."""
    every = read_rows(path)[1]
    rows = distinct(every)
    count = {}
    for r in every:
        key = r.replace("-", "")
        count[key] = count.get(key, 0) + 1
    weight = [count[r.replace("-", "")] for r in rows]
    unit = refine = affine = 0
    for a in range(len(rows)):
        for b in range(a + 1, len(rows)):
            w = weight[a] * weight[b]
            d, m, g, o = pair(rows[a], rows[b])
            unit += w * d
            refine += w * (d + o)
            affine += w * (4 * m + 2 * g + 6 * o)
    return unit, refine, affine


def one(task):
    ours, theirs = task[1], task[2]
    if not (os.path.exists(ours) and os.path.exists(theirs)):
        return task[0], None, None
    try:
        return task[0], costs(ours), costs(theirs)
    except Exception:
        return task[0], None, None


def metric(ve, method, rid):
    p = os.path.join(ve, "results", method, rid + ".json")
    if not os.path.exists(p):
        return None
    a = json.load(open(p)).get("alignment", {})
    u, f = a.get("all_cost_over_opt"), a.get("affine_all_cost_over_opt")
    return None if u is None or f is None else (u, f)


def main():
    args = sys.argv[1:]
    jobs = 8
    if "--jobs" in args:
        k = args.index("--jobs")
        jobs = int(args[k + 1])
        del args[k:k + 2]
    ve, method, other_dir, other = args
    table = list(csv.DictReader(open(os.path.join(ve, "regions", "regions.tsv")), delimiter="\t"))
    info = {r["region_id"]: r for r in table}
    tasks = [(rid, os.path.join(ve, "work", "plegma", method, rid + ".msa.fa"), os.path.join(other_dir, rid + ".msa.fa"))
             for rid in info]
    with Pool(jobs) as pool:
        got = pool.map(one, tasks, chunksize=1)
    classes = collections.defaultdict(collections.Counter)
    lost = []
    for rid, mine, theirs in got:
        a, b = metric(ve, method, rid), metric(ve, other, rid)
        if mine is None or theirs is None or a is None or b is None:
            continue
        st = info[rid]["stratum"]
        for k, name in ((0, "unit"), (1, "affine")):
            if b[k] < a[k] - 0.002:
                cls = "search" if theirs[1] < mine[1] else "objective"
                classes[st]["%s lost: %s" % (name, cls)] += 1
                if k == 0:
                    lost.append((a[0] - b[0], rid, st, cls, a, b, mine, theirs))
            elif a[k] < b[k] - 0.002:
                classes[st]["%s won" % name] += 1
            else:
                classes[st]["%s tied" % name] += 1
    print("per stratum (unit and affine metrics; a loss is `search` when the other MSA is lower on plegma's own"
          " objective):")
    for st in sorted(classes):
        print("  %-21s %s" % (st, dict(sorted(classes[st].items()))))
    lost.sort(reverse=True)
    print("\nlost on the unit metric, by most (ours / other: unit, affine metrics; objective; unit SP; affine SP):")
    for d, rid, st, cls, a, b, mine, theirs in lost[:25]:
        print("  %-9s %-21s %-9s %.3f/%.3f vs %.3f/%.3f  obj %d vs %d (%+.2f%%)  unit %d vs %d  affine %d vs %d" % (
            rid, st, cls, a[0], a[1], b[0], b[1], mine[1], theirs[1], 100.0 * (theirs[1] - mine[1]) / mine[1],
            mine[0], theirs[0], mine[2], theirs[2]))


if __name__ == "__main__":
    main()
