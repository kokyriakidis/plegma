#!/usr/bin/env python3
"""One iteration of plegma on vntr-eval: candidates, scores and the comparison, as fast as the machine allows.

usage: python3 -I experiments/vntreval/iterate.py <vntr-eval> <plegma> <label> [--flags "..."] [--jobs n]
           [--quick | --all | --regions a,b,...] [--against m,m,...] [--truth]

- The method is named `<label>_<first 10 hex of the binary's sha256 and the flags'>`: a binary and flags already
  scored are not run again, and each version's results stay side by side under results/.
- Candidates: `plegma <flags> -r 1 <region>/hap32.fa` (`-L -I 20,64,1` by default), then vntr-eval's
  msa_graph.py, every region at once over `--jobs` workers (8), longest regions first.
- Scores: vntr-eval's evaluate_all.py, without the truth stage unless `--truth` (the pairwise optima are cached
  under vntr-eval's work/cache).
- Regions: `--quick` (24: where plegma loses most to the unit-aware realigner, the hotspots traced so far, six
  controls), the development set (every hotspot and 16 controls) by default, `--all` for the 149, or `--regions`.
- Prints compare.py's tables, unit and affine, against abPOA (`abpoa_longest`), pgrealign's medoid star
  (`pgrealign_star`, star_candidates.py), the unit-aware realigner and the methods named in `--against`, and the
  wall time of each step.
"""
import csv
import hashlib
import os
import subprocess
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
QUICK = ["L011463", "L000308", "L002143", "L003872", "L010060", "L007172", "L011224", "L008902", "L013064",
         "L008216", "L015223", "L007161", "L001909", "L009658", "L015904", "L008877", "L012256", "L007081",
         "L009283", "L007065", "L001899", "L003883", "L013399", "L004620"]
"""A minute or two: where plegma loses most to the unit-aware realigner, the hotspots traced so far, six controls."""
DEV_CONTROLS = ["L009283", "L003883", "L003872", "L007065", "L001899", "L000381", "L015223", "L012313",
                "L001491", "L008902", "L000407", "L003672", "L015285", "L014426", "L013399", "L004620"]


def option(args, name, default=None, flag=False):
    if name in args:
        k = args.index(name)
        if flag:
            del args[k]
            return True
        value = args[k + 1]
        del args[k:k + 2]
        return value
    return False if flag else default


def candidate(task):
    vntr_eval, plegma, flags, method, rid = task
    msa_dir = os.path.join(vntr_eval, "work", "plegma", method)
    gfa = os.path.join(vntr_eval, "candidates", method, rid + ".gfa")
    if os.path.exists(gfa) and os.path.getsize(gfa) > 0:
        return rid, 0.0, "cached"
    hap = os.path.join(vntr_eval, "regions", rid, "hap32.fa")
    msa = os.path.join(msa_dir, rid + ".msa.fa")
    start = time.time()
    with open(msa, "w") as out:
        r = subprocess.run([plegma] + flags + ["-r", "1", hap], stdout=out, stderr=subprocess.PIPE, timeout=7200)
    seconds = time.time() - start
    if r.returncode != 0:
        return rid, seconds, "plegma failed: " + r.stderr.decode()[-200:]
    g = subprocess.run(["python3", "-I", os.path.join(vntr_eval, "tools", "msa_graph.py"), msa, hap, gfa, "--quiet"],
                       capture_output=True, text=True)
    return rid, seconds, "ok" if g.returncode == 0 else "graph failed: " + g.stderr[-200:]


def main():
    args = sys.argv[1:]
    flags = option(args, "--flags", "-L -I 20,64,1").split()
    jobs = int(option(args, "--jobs", "8"))
    every = option(args, "--all", flag=True)
    quick = option(args, "--quick", flag=True)
    chosen = option(args, "--regions")
    against = [m for m in option(args, "--against", "").split(",") if m]
    truth = option(args, "--truth", flag=True)
    vntr_eval, plegma, label = args
    h = hashlib.sha256(open(plegma, "rb").read())
    h.update(" ".join(flags).encode())
    method = "%s_%s" % (label, h.hexdigest()[:10])
    os.makedirs(os.path.join(vntr_eval, "candidates", method), exist_ok=True)
    os.makedirs(os.path.join(vntr_eval, "work", "plegma", method), exist_ok=True)
    with open(os.path.join(vntr_eval, "regions", "regions.tsv")) as f:
        table = list(csv.DictReader(f, delimiter="\t"))
    if chosen:
        regions = chosen.split(",")
    elif quick:
        regions = list(QUICK)
    elif every:
        regions = [r["region_id"] for r in table]
    else:
        regions = [r["region_id"] for r in table if r["stratum"].startswith("hotspot")] + DEV_CONTROLS
    size = {r["region_id"]: int(r["span_end"]) - int(r["span_start"]) for r in table}
    regions.sort(key=lambda rid: -size.get(rid, 0))
    print("method %s: %d regions, plegma %s" % (method, len(regions), " ".join(flags)), flush=True)
    t0 = time.time()
    with Pool(jobs) as pool:
        results = pool.map(candidate, [(vntr_eval, plegma, flags, method, rid) for rid in regions], chunksize=1)
    failed = [(rid, msg) for rid, _, msg in results if msg not in ("ok", "cached")]
    ran = [s for _, s, msg in results if msg == "ok"]
    print("candidates: %.0f s wall, %d run (plegma %.0f s summed, slowest %.0f s), %d cached, %d failed" % (
        time.time() - t0, len(ran), sum(ran), max(ran, default=0), sum(m == "cached" for _, _, m in results),
        len(failed)), flush=True)
    for rid, msg in failed:
        print("  %s: %s" % (rid, msg))
    t1 = time.time()
    cmd = ["python3", "-I", os.path.join(vntr_eval, "tools", "evaluate_all.py"), method, "--jobs", str(jobs),
           "--threads", "1", "--regions", ",".join(regions)]
    if not truth:
        cmd += ["--", "--no-truvari"]
    subprocess.run(cmd, cwd=vntr_eval, capture_output=True, text=True)
    print("scores: %.0f s wall" % (time.time() - t1), flush=True)
    compare = os.path.join(HERE, "compare.py")
    for metric in ("alignment.all_cost_over_opt", "alignment.affine_all_cost_over_opt"):
        print("\n## %s" % metric)
        out = subprocess.run(["python3", "-I", compare, vntr_eval, "abpoa_longest", method, "pgrealign_star",
                              "unit_aware"] + against +
                             ["--metric", metric, "--worst", "0"], capture_output=True, text=True).stdout
        print("\n".join(l for l in out.splitlines() if not l.startswith("Where")))
    print("method: %s" % method)


if __name__ == "__main__":
    main()
