#!/usr/bin/env python3
"""The anchor audit (build/anchors) over a panel of vntr-eval regions: each plan's off anchors and letters, by region
and summed.

usage: python3 -I experiments/anchors/panel.py <vntr-eval> <anchors binary> <placement binary> <out dir>
           [--regions a,b,...] [--pairs n] [--jobs n] -- <plan> [<plan> ...]

Each region's period is `placement --tandem`'s (its copy layout's, else 0). Writes <out dir>/<region>.tsv (the
anchors) and prints, per plan, the regions' anchored letters and the share off, region by region and over all.
"""
import os
import re
import subprocess
import sys
from multiprocessing import Pool

QUICK = ["L011463", "L000308", "L002143", "L003872", "L010060", "L007172", "L011224", "L008902", "L013064",
         "L008216", "L015223", "L007161", "L001909", "L009658", "L015904", "L008877", "L012256", "L007081",
         "L009283", "L007065", "L001899", "L003883", "L013399", "L004620"]


def one(task):
    ve, anchors, placement, out, rid, pairs, plans = task
    hap = os.path.join(ve, "regions", rid, "hap32.fa")
    tandem = subprocess.run([placement, "--tandem", hap], capture_output=True, text=True).stdout
    m = re.search(r"^layout period (\d+)", tandem, re.M)
    period = int(m.group(1)) if m else 0
    run = subprocess.run([anchors, hap, str(period), str(pairs), os.path.join(out, rid + ".tsv")] + plans,
                         capture_output=True, text=True)
    totals = {}
    for line in run.stdout.splitlines():
        m = re.match(r"^(\S+) \| anchors (\d+) off (\d+) \| letters (\d+) off (\d+).*forced excess (-?\d+) of optima"
                     r" (\d+)", line)
        if m:
            totals[m.group(1)] = tuple(int(m.group(k)) for k in range(2, 8))
    return rid, period, totals, run.stderr[-300:] if run.returncode else ""


def main():
    args = sys.argv[1:]
    k = args.index("--")
    plans = args[k + 1:]
    args = args[:k]

    def option(name, default):
        if name in args:
            i = args.index(name)
            v = args[i + 1]
            del args[i:i + 2]
            return v
        return default

    regions = option("--regions", ",".join(QUICK)).split(",")
    pairs = int(option("--pairs", "30"))
    jobs = int(option("--jobs", "8"))
    ve, anchors, placement, out = args
    os.makedirs(out, exist_ok=True)
    with Pool(jobs) as pool:
        got = pool.map(one, [(ve, anchors, placement, out, rid, pairs, plans) for rid in regions], chunksize=1)
    sums = {p: [0, 0, 0, 0, 0, 0] for p in plans}
    print("%-9s %6s  %s" % ("region", "period", "  ".join("%-22s" % p[:22] for p in plans)))
    for rid, period, totals, err in got:
        cells = []
        for p in plans:
            t = totals.get(p)
            if t is None:
                cells.append("%-22s" % "-")
                continue
            for i in range(6):
                sums[p][i] += t[i]
            cells.append("%-22s" % ("%5.1f%% +%.2f%%" % (100.0 * t[3] / max(1, t[2]), 100.0 * t[4] / max(1, t[5]))))
        print("%-9s %6d  %s%s" % (rid, period, "  ".join(cells), ("  ERROR " + err) if err else ""))
    print("\ncells: letters off an optimal alignment, then the excess the anchors force over the pairs' optima")
    print("\nover all (anchors, off anchors, letters, off letters, off share, forced excess / optima):")
    for p in plans:
        s = sums[p]
        print("  %-48s %7d %7d %9d %9d %6.2f%%  +%.3f%%" % (p, s[0], s[1], s[2], s[3], 100.0 * s[3] / max(1, s[2]),
                                                         100.0 * s[4] / max(1, s[5])))


if __name__ == "__main__":
    main()
