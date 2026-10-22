#!/usr/bin/env python3
"""Runs build/seeding over a pair set in parallel chunks and summarises each strategy.

usage: python3 -I experiments/seeding/run.py <pairs.fa> <truth.tsv or -> <out dir> <jobs> <strategy> ... [--meta file]

With `--meta` (real_pairs.py's), results are also broken down by stratum and by repeat period.

Writes <out dir>/pairs.tsv (every pair under every strategy) and prints, per strategy: pairs, mean F1 against the
truth and the share aligned exactly as the truth, mean F1 against the optimum, cost over optimum (pooled) and
the share at the optimum, anchors per pair and the share lying on the truth, the largest window (median) and
time; then the same F1 and exactness by repeat period, copy variation and copy-number difference.
"""
import collections
import csv
import os
import statistics
import subprocess
import sys
from multiprocessing import Pool

BINARY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "build", "seeding")


def chunk_run(task):
    fasta, truth, out, strategies = task
    subprocess.run([BINARY, fasta, truth, out] + strategies, check=True)
    return out


def main():
    pairs, truth, out_dir, jobs = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    strategies = sys.argv[5:]
    meta_file = None
    if "--meta" in strategies:
        meta_file = strategies[strategies.index("--meta") + 1]
        strategies = strategies[:strategies.index("--meta")]
    os.makedirs(out_dir, exist_ok=True)
    records = []
    with open(pairs) as f:
        for line in f:
            if line.startswith(">"):
                records.append([line[1:].strip(), ""])
            else:
                records[-1][1] += line.strip()
    meta_lines = open(truth).read().splitlines() if truth != "-" else []
    header = meta_lines[0] if meta_lines else ""
    n = len(records) // 2
    per = max(1, (n + jobs - 1) // jobs)
    tasks = []
    for c in range(0, n, per):
        fa = os.path.join(out_dir, "chunk%d.fa" % c)
        with open(fa, "w") as f:
            for k in range(c, min(n, c + per)):
                for name, seq in records[2 * k: 2 * k + 2]:
                    f.write(">%s\n%s\n" % (name, seq))
        tr = "-"
        if meta_lines:
            tr = os.path.join(out_dir, "chunk%d.truth.tsv" % c)
            with open(tr, "w") as f:
                f.write(header + "\n" + "\n".join(meta_lines[1 + c: 1 + min(n, c + per)]) + "\n")
        tasks.append((fa, tr, os.path.join(out_dir, "chunk%d.out.tsv" % c), strategies))
    with Pool(jobs) as pool:
        outs = pool.map(chunk_run, tasks)
    rows = []
    for o in outs:
        rows.extend(csv.DictReader(open(o), delimiter="\t"))
    with open(os.path.join(out_dir, "pairs.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    meta = {}
    for line in meta_lines[1:]:
        fields = line.split("\t")
        meta[fields[0]] = fields
    has_truth = bool(meta)
    by = collections.defaultdict(list)
    for r in rows:
        by[r["strategy"]].append(r)
    print("| strategy | pairs | F1 truth | exact | F1 optimum | cost/opt | at optimum | anchors | on truth | window | ms |")
    print("| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |")
    for s in strategies:
        g = by[s]
        f1 = [float(r["f1_truth"]) for r in g]
        exact = sum(x >= 0.999999 for x in f1) / len(g)
        anchors = sum(int(r["anchors"]) for r in g)
        true_anchors = sum(int(r["anchors_true"]) for r in g)
        print("| %s | %d | %s | %s | %.4f | %.4f | %.2f | %.1f | %s | %d | %.1f |" % (
            s, len(g),
            "%.4f" % statistics.mean(f1) if has_truth else "-",
            "%.2f" % exact if has_truth else "-",
            statistics.mean(float(r["f1_optimum"]) for r in g),
            sum(int(r["cost"]) for r in g) / max(1, sum(int(r["optimum"]) for r in g)),
            sum(r["cost"] == r["optimum"] for r in g) / len(g),
            anchors / len(g),
            "%.3f" % (true_anchors / anchors) if has_truth and anchors else "-",
            statistics.median(int(r["largest_window"]) for r in g),
            sum(float(r["ms"]) for r in g)))
    if meta_file:
        info = {}
        for line in open(meta_file).read().splitlines()[1:]:
            fields = line.split("\t")
            info[fields[0]] = fields
        def period_class(m):
            try:
                p = float(m[3])
            except ValueError:
                return "?"
            return "<=6" if p <= 6 else ("<=30" if p <= 30 else ">30")
        for label, key in [("stratum", lambda m: m[2]), ("period", period_class), ("pair kind", lambda m: m[5])]:
            groups = sorted({key(m) for m in info.values()})
            print("\n| %s | %s |" % (label, " | ".join("%s (F1 opt / cost-opt / at opt)" % g for g in groups)))
            print("| :-- |" + " --: |" * len(groups))
            for s in strategies:
                cells = []
                for gname in groups:
                    g = [r for r in by[s] if key(info[r["pair"][:-2]]) == gname]
                    if not g:
                        cells.append("-")
                        continue
                    cells.append("%.4f / %.4f / %.2f (%d)" % (
                        statistics.mean(float(r["f1_optimum"]) for r in g),
                        sum(int(r["cost"]) for r in g) / max(1, sum(int(r["optimum"]) for r in g)),
                        sum(r["cost"] == r["optimum"] for r in g) / len(g), len(g)))
                print("| %s | %s |" % (s, " | ".join(cells)))
    if has_truth:
        optimum_f1 = statistics.mean(float(r["optimum_f1_truth"]) for r in by[strategies[0]])
        exact_opt = sum(float(r["optimum_f1_truth"]) >= 0.999999 for r in by[strategies[0]]) / len(by[strategies[0]])
        print("\nThe optimal alignment itself: F1 truth %.4f, exact %.2f.\n" % (optimum_f1, exact_opt))
        for label, key in [("period", lambda m: "<=3" if int(m[1]) <= 3 else ("<=18" if int(m[1]) <= 18 else ">18")),
                           ("copy variation", lambda m: m[3]),
                           ("copy-number difference", lambda m: "0" if int(m[5]) == 0 else ("<=5" if abs(int(m[5])) <= 5 else ">5"))]:
            groups = sorted({key(m) for m in meta.values()})
            print("| %s | %s |" % (label, " | ".join("%s (F1 / exact)" % g for g in groups)))
            print("| :-- |" + " --: |" * len(groups))
            for s in ["optimum"] + strategies:
                cells = []
                for gname in groups:
                    if s == "optimum":
                        g = [r for r in by[strategies[0]] if key(meta[r["pair"][:-2]]) == gname]
                        vals = [float(r["optimum_f1_truth"]) for r in g]
                    else:
                        g = [r for r in by[s] if key(meta[r["pair"][:-2]]) == gname]
                        vals = [float(r["f1_truth"]) for r in g]
                    cells.append("%.4f / %.2f" % (statistics.mean(vals), sum(v >= 0.999999 for v in vals) / len(vals)) if vals else "-")
                print("| %s | %s |" % (s, " | ".join(cells)))
            print()


if __name__ == "__main__":
    main()
