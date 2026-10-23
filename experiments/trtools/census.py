#!/usr/bin/env python3
"""longdust and TRF over every vntr-eval region: does a tandem period explain each region's low-complexity array?

usage: python3 -I experiments/trtools/census.py <vntr-eval> <out dir> [--jobs n] [--compare a,b]

For every region, its distinct haplotypes go through longdust (-k8 -w20000) and TRF (2 7 7 80 10 50 2000); per
haplotype, longdust's longest interval is the array, and TRF's best period inside it is the smallest period whose repeats
cover within 10% as much of it as the period that covers most (TRF reports multiples of a period too); TRF's output
is kept in the out dir and reused. The region's class:
- `tandem`: in most haplotypes one period's repeats cover at least 80% of the array, the same period in all of them;
- `drifting`: they cover at least 80%, but the period differs between haplotypes;
- `weak`: an array, but no period covers 80% of it (low complexity without clear tandem structure);
- `none`: no longdust interval of 200 letters or more.
Prints each region's class, the annotated period and TRF's, the coverage, and (with `--compare a,b`) the two methods'
all_cost_over_opt; then the classes per stratum, and the methods' medians per class. The tools come from
$TRTOOLS/bin.
"""
import collections
import csv
import json
import os
import statistics
import subprocess
import sys
from multiprocessing import Pool

B = os.path.join(os.environ.get("TRTOOLS", os.path.expanduser("~/micromamba/envs/trtools")), "bin")


def distinct_fasta(src, dst):
    names, seqs = [], []
    for line in open(src):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0].replace("#", "_"))
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    seen = set()
    with open(dst, "w") as f:
        for n, s in zip(names, ("".join(x) for x in seqs)):
            if s not in seen:
                seen.add(s)
                f.write(">%s\n%s\n" % (n, s))
    return len(seen)


def one(task):
    ve, out, rid = task
    d = os.path.join(out, rid)
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, "hap.fa")
    distinct_fasta(os.path.join(ve, "regions", rid, "hap32.fa"), fa)
    ld = subprocess.run([os.path.join(B, "longdust"), "-k8", "-w20000", fa], capture_output=True, text=True).stdout
    if not os.path.exists(os.path.join(d, "trf.dat")):
        with open(os.path.join(d, "trf.part"), "w") as f:
            subprocess.run([os.path.join(B, "trf"), "hap.fa", "2", "7", "7", "80", "10", "50", "2000", "-h", "-ngs"],
                           cwd=d, stdout=f, stderr=subprocess.DEVNULL)
        os.replace(os.path.join(d, "trf.part"), os.path.join(d, "trf.dat"))
    arrays = {}
    for line in ld.splitlines():
        c = line.split("\t")
        s, e = int(c[1]), int(c[2])
        if e - s >= 200 and (c[0] not in arrays or e - s > arrays[c[0]][1] - arrays[c[0]][0]):
            arrays[c[0]] = (s, e)
    repeats = collections.defaultdict(list)
    name = None
    for line in open(os.path.join(d, "trf.dat")):
        if line.startswith("@"):
            name = line[1:].strip()
        elif name and line.strip():
            c = line.split()
            repeats[name].append((int(c[0]) - 1, int(c[1]), int(c[2])))
    per_hap = []
    for h, (s, e) in arrays.items():
        # Each period's coverage of the array (its repeats' union inside it).
        cover = {}
        for p in {r[2] for r in repeats[h]}:
            marks = bytearray(e - s)
            for a, b, q in repeats[h]:
                if q == p:
                    for x in range(max(a, s), min(b, e)):
                        marks[x - s] = 1
            cover[p] = sum(marks) / (e - s)
        # TRF also reports a period's multiples, which cover as much: the smallest period within 10% of the best.
        top = max(cover.values()) if cover else 0.0
        best = min(((p, c) for p, c in cover.items() if c >= 0.9 * top), default=(0, 0.0))
        per_hap.append((h, e - s, best[0], best[1]))
    return rid, per_hap


def classify(per_hap):
    if not per_hap:
        return "none", 0, 0.0
    good = [h for h in per_hap if h[3] >= 0.8]
    periods = collections.Counter(h[2] for h in good)
    median_cover = statistics.median(h[3] for h in per_hap)
    if len(good) * 2 > len(per_hap):
        top, count = periods.most_common(1)[0]
        # Periods within 5% of the commonest count as the same.
        same = sum(n for p, n in periods.items() if abs(p - top) <= 0.05 * top)
        return ("tandem" if same == len(good) else "drifting"), top, median_cover
    top = collections.Counter(h[2] for h in per_hap).most_common(1)[0][0]
    return "weak", top, median_cover


def main():
    args = sys.argv[1:]
    jobs = 8
    compare = []
    if "--jobs" in args:
        k = args.index("--jobs")
        jobs = int(args[k + 1])
        del args[k:k + 2]
    if "--compare" in args:
        k = args.index("--compare")
        compare = args[k + 1].split(",")
        del args[k:k + 2]
    ve, out = args
    table = list(csv.DictReader(open(os.path.join(ve, "regions", "regions.tsv")), delimiter="\t"))
    info = {r["region_id"]: r for r in table}
    with Pool(jobs) as pool:
        results = dict(pool.map(one, [(ve, out, r["region_id"]) for r in table], chunksize=1))

    def score(m, rid):
        p = os.path.join(ve, "results", m, rid + ".json")
        if not os.path.exists(p):
            return None
        return json.load(open(p)).get("alignment", {}).get("all_cost_over_opt")

    by_stratum = collections.defaultdict(collections.Counter)
    by_class = collections.defaultdict(lambda: collections.defaultdict(list))
    print("%-9s %-21s %-9s %8s %6s %6s  %s" % ("region", "stratum", "class", "annotated", "TRF", "cover",
                                             "  ".join(compare)))
    for rid in sorted(results, key=lambda r: (info[r]["stratum"], r)):
        cls, period, cover = classify(results[rid])
        st = info[rid]["stratum"]
        by_stratum[st][cls] += 1
        scores = [score(m, rid) for m in compare]
        for m, v in zip(compare, scores):
            if v is not None:
                by_class[cls][m].append(v)
        print("%-9s %-21s %-9s %8s %6s %6.2f  %s" % (rid, st, cls, info[rid].get("period", ""), period, cover,
                                                    "  ".join("%.3f" % v if v else "-" for v in scores)))
    print("\nclasses per stratum:")
    for st, c in sorted(by_stratum.items()):
        print("  %-21s %s" % (st, dict(c)))
    if compare:
        print("\nmedian all_cost_over_opt per class:")
        for cls, ms in sorted(by_class.items()):
            print("  %-9s %3d regions  %s" % (cls, len(next(iter(ms.values()))), "  ".join(
                "%s %.4f" % (m, statistics.median(v)) for m, v in ms.items())))


if __name__ == "__main__":
    main()
