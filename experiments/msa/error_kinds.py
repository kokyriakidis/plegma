#!/usr/bin/env python3
"""Error kinds of every consensus accuracy.py left, per configuration and stratum.

usage: python3 -I experiments/msa/error_kinds.py <karan> <set dir> <jobs> name ...

For every `<set>.fa` with a `<set>.truth.fa` in `set dir` and every configuration `name` whose consensus
`<set>.fa.<name>.cons.fa` exists, sorts its errors (errors.py) and prints, per configuration, the totals by kind
over all sets, and by stratum (the set name's second field: vntr_<stratum>_...).
"""
import collections
import glob
import os
import subprocess
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
KINDS = ["substitution", "homopolymer", "tandem", "other"]


def classify(task):
    karan, truth, cons, name, stratum = task
    out = subprocess.run(["python3", "-I", os.path.join(HERE, "errors.py"), karan, truth, cons],
                         capture_output=True, text=True).stdout.split()
    return name, stratum, [int(x) for x in out[1:5]] if len(out) == 5 else None


def main():
    karan, folder, jobs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    names = sys.argv[4:]
    tasks = []
    for fasta in sorted(glob.glob(os.path.join(folder, "*.fa"))):
        if fasta.endswith(".truth.fa") or ".cons." in fasta:
            continue
        truth = fasta[:-3] + ".truth.fa"
        base = os.path.basename(fasta)
        stratum = base.split("_")[1] + "_" + base.split("_")[2] if base.startswith("vntr_") else "all"
        for name in names:
            cons = fasta + "." + name + ".cons.fa"
            if os.path.exists(truth) and os.path.exists(cons):
                tasks.append((karan, truth, cons, name, stratum))
    with Pool(jobs) as pool:
        results = pool.map(classify, tasks)
    total = collections.defaultdict(lambda: [0] * 4)
    by_stratum = collections.defaultdict(lambda: [0] * 4)
    failed = collections.Counter()
    for name, stratum, counts in results:
        if counts is None:
            failed[name] += 1
            continue
        for k in range(4):
            total[name][k] += counts[k]
            by_stratum[(name, stratum)][k] += counts[k]
    print("| config | " + " | ".join(KINDS) + " | non-homopolymer | failed |")
    print("| :-- |" + " --: |" * 6)
    for name in names:
        t = total[name]
        print("| %s | %s | %d | %d |" % (name, " | ".join(str(x) for x in t), t[0] + t[2] + t[3], failed[name]))
    strata = sorted({s for _, s in by_stratum})
    print("\n| config | " + " | ".join("%s (subst / tandem / other)" % s for s in strata) + " |")
    print("| :-- |" + " --: |" * len(strata))
    for name in names:
        cells = []
        for s in strata:
            c = by_stratum.get((name, s), [0] * 4)
            cells.append("%d / %d / %d" % (c[0], c[2], c[3]))
        print("| %s | %s |" % (name, " | ".join(cells)))


if __name__ == "__main__":
    main()
