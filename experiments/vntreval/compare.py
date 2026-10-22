#!/usr/bin/env python3
"""Methods against a rival on vntr-eval's results, region by region.

usage: python3 -I experiments/vntreval/compare.py <vntr-eval> <rival> <method> ... [--metric name] [--worst n]

For each stratum, on the regions where the rival and every method have a valid result: the median of the metric
(`alignment.all_cost_over_opt` by default: the graph-implied pairwise cost over the optimum, lower better) for
each, and for each method how many regions it is better, worse or equal on than the rival, and its median
difference. Then the `n` regions (10 by default) where the first method loses most to the rival.
"""
import collections
import glob
import json
import os
import statistics
import sys


def value(vntr_eval, method, region, metric):
    path = os.path.join(vntr_eval, "results", method, region + ".json")
    if not os.path.exists(path):
        return None, None
    data = json.load(open(path))
    if data.get("status") != "ok":
        return data.get("stratum"), None
    node = data
    for part in metric.split("."):
        node = node.get(part) if isinstance(node, dict) else None
    return data.get("stratum"), node


def main():
    args = sys.argv[1:]
    metric = args.pop(args.index("--metric") + 1) if "--metric" in args else "alignment.all_cost_over_opt"
    if "--metric" in args:
        args.remove("--metric")
    worst = int(args.pop(args.index("--worst") + 1)) if "--worst" in args else 10
    if "--worst" in args:
        args.remove("--worst")
    vntr_eval, rival, methods = args[0], args[1], args[2:]
    regions = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(vntr_eval, "results", rival, "*.json")))
    rows = collections.defaultdict(list)
    for region in regions:
        stratum, base = value(vntr_eval, rival, region, metric)
        values = [value(vntr_eval, m, region, metric)[1] for m in methods]
        if base is None or any(v is None for v in values):
            continue
        rows[stratum].append((region, base, values))
    print("| stratum | regions | %s | %s |" % (rival, " | ".join("%s (better / worse / equal, median diff)" % m
                                                              for m in methods)))
    print("| :-- | --: |" + " --: |" * (len(methods) + 1))
    for stratum in sorted(rows):
        r = rows[stratum]
        cells = []
        for k in range(len(methods)):
            d = [x[2][k] - x[1] for x in r]
            cells.append("%.4f (%d / %d / %d, %+.4f)" % (statistics.median(x[2][k] for x in r),
                                                          sum(v < -1e-9 for v in d), sum(v > 1e-9 for v in d),
                                                          sum(abs(v) <= 1e-9 for v in d), statistics.median(d)))
        print("| %s | %d | %.4f | %s |" % (stratum, len(r), statistics.median(x[1] for x in r), " | ".join(cells)))
    every = sorted((x for r in rows.values() for x in r), key=lambda x: x[1] - x[2][0])[:worst]
    print("\nWhere %s loses most to %s:" % (methods[0], rival))
    for region, base, values in every:
        print("  %s %s %.3f %s" % (region, rival, base, " ".join("%s %.3f" % (m, v) for m, v in zip(methods, values))))


if __name__ == "__main__":
    main()
