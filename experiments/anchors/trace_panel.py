#!/usr/bin/env python3
"""plegma's anchors (split_anchors) traced against the optimum over the brute-force panel's pairs
(build/anchors --brute ... --anchors): every anchor off the optimum, with its stage, region, period and how it is off.

usage: python3 -I experiments/anchors/trace_panel.py <vntr-eval> <anchors binary> <placement binary> <panel dir>
           <out.tsv> [--prices x,o,e] [--jobs n] [--ambiguity share]

The pairs are those <panel dir> (brute_panel.py's) holds, <region>.<a>_<b>.txt.
"""
import os
import re
import subprocess
import sys
from multiprocessing import Pool


def one(task):
    ve, anchors, placement, rid, a, b, prices, ambiguity = task
    hap = os.path.join(ve, "regions", rid, "hap32.fa")
    tandem = subprocess.run([placement, "--tandem", hap], capture_output=True, text=True).stdout
    m = re.search(r"^layout period (\d+)", tandem, re.M)
    period = int(m.group(1)) if m else 0
    out = subprocess.run([anchors, "--brute", hap, str(a), str(b), str(max(period, 1)), "--prices", prices,
                          "--anchors", "--ambiguity", ambiguity], capture_output=True, text=True).stdout
    return rid, a, b, period, out


def main():
    args = sys.argv[1:]
    prices, jobs, ambiguity = "1,0,1", 9, "0"
    if "--ambiguity" in args:
        ambiguity = args[args.index("--ambiguity") + 1]
        del args[args.index("--ambiguity"):args.index("--ambiguity") + 2]
    if "--prices" in args:
        prices = args[args.index("--prices") + 1]
        del args[args.index("--prices"):args.index("--prices") + 2]
    if "--jobs" in args:
        jobs = int(args[args.index("--jobs") + 1])
        del args[args.index("--jobs"):args.index("--jobs") + 2]
    ve, anchors, placement, panel, out = args
    tasks = []
    for f in sorted(os.listdir(panel)):
        m = re.match(r"^(\S+)\.(\d+)_(\d+)\.txt$", f)
        if m:
            tasks.append((ve, anchors, placement, m.group(1), int(m.group(2)), int(m.group(3)), prices, ambiguity))
    keys = ["stage", "a", "b", "len", "excess", "right", "lead", "trail", "trimmed", "diag"]
    with Pool(jobs) as pool, open(out, "w") as f:
        f.write("region\tpa\tpb\tperiod\t" + "\t".join(keys) + "\toptimal\tlengths\n")
        joint = open(out + ".joint", "w")
        for rid, a, b, period, text in pool.imap_unordered(one, tasks):
            for line in text.splitlines():
                if line.startswith("joint "):
                    w = line.split()
                    joint.write("%s\t%d\t%d\t%s\t%s\n" % (rid, a, b, w[1], w[3]))
                if not line.startswith("anchor "):
                    continue
                w = line.split()
                rec = {"stage": w[1], "a": w[3], "b": w[5], "len": w[7], "excess": w[9], "right": w[11],
                       "lead": w[13], "trail": w[15], "trimmed": w[17], "diag": w[19]}
                optimal = " ".join(w[21:27])
                lengths = " ".join(w[28:30])
                f.write("%s\t%d\t%d\t%d\t%s\t%s\t%s\n" % (rid, a, b, period, "\t".join(rec[k] for k in keys), optimal,
                                                          lengths))


if __name__ == "__main__":
    main()
