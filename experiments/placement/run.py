#!/usr/bin/env python3
"""The placement harness over many regions: every row of each region's MSA placed afresh by every strategy.

usage: python3 -I experiments/placement/run.py <placement binary> <msa dir> <out dir> <jobs> <regions> <strategy> ...

<msa dir> holds `<region>.msa.fa` (an MSA of the region's hap32 sequences, e.g. a vntr-eval candidate's); <regions>
a comma list. Writes <out dir>/<region>.tsv per region and prints, per strategy: the placed rows' cost over their
floor (pooled: summed placed over summed floor), the present rows' over the same floor, the share of rows the
placement lowered and raised, and the median milliseconds a placement took.
"""
import csv
import os
import statistics
import subprocess
import sys
from multiprocessing import Pool


def run(task):
    binary, msa_dir, out_dir, rid, strategies = task
    out = os.path.join(out_dir, rid + ".tsv")
    if not os.path.exists(out):
        subprocess.run([binary, os.path.join(msa_dir, rid + ".msa.fa"), rid, out] + strategies, check=True,
                       capture_output=True)
    return out


def main():
    binary, msa_dir, out_dir, jobs, regions = sys.argv[1:6]
    strategies = sys.argv[6:]
    os.makedirs(out_dir, exist_ok=True)
    with Pool(int(jobs)) as pool:
        files = pool.map(run, [(binary, msa_dir, out_dir, rid, strategies) for rid in regions.split(",")], chunksize=1)
    rows = []
    for f in files:
        rows += list(csv.DictReader(open(f), delimiter="\t"))
    print("| strategy | placed / floor | present / floor | rows lowered | rows raised | median ms |")
    print("| :-- | --: | --: | --: | --: | --: |")
    for s in strategies:
        rs = [r for r in rows if r["strategy"] == s]
        placed = sum(int(r["placed"]) for r in rs)
        present = sum(int(r["present"]) for r in rs)
        floor = sum(int(r["floor"]) for r in rs)
        lower = sum(int(r["placed"]) < int(r["present"]) for r in rs)
        higher = sum(int(r["placed"]) > int(r["present"]) for r in rs)
        print("| %s | %.4f | %.4f | %.2f | %.2f | %.0f |" % (s, placed / floor, present / floor, lower / len(rs),
                                                         higher / len(rs), statistics.median(float(r["ms"]) for r in rs)))


if __name__ == "__main__":
    main()
