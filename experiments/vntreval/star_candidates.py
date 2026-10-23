#!/usr/bin/env python3
"""pgrealign's medoid star on every vntr-eval region package, as candidates/<method>/<id>.gfa for evaluate_all.py.

usage: python3 -I experiments/vntreval/star_candidates.py <vntr-eval> <pgrealign repository> <abpoa> [method]
           [--jobs n]

Each region's hap32 sequences are its alleles, each weighted by how many of the 34 spell it (as pgrealign weights
an allele by its walks), and aligned by pgrealign's own star_msa with realign's defaults: the medoid by k-mer
content as the centre, abPOA (-M 8 -X 16 -O 16,96 -E 8,1, pgrealign's scores) for each allele against it, the
distinct insertions of one slot aligned to each other, `max_cells` as realign's; no homology guard (every region is
aligned). The MSA goes through vntr-eval's msa_graph.py, as every method's does. The method is `pgrealign_star` by
default; a region whose graph exists is skipped.
"""
import collections
import os
import subprocess
import sys
import time
from multiprocessing import Pool

PGREALIGN = ABPOA = None


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return names, ["".join(s) for s in seqs]


def one(task):
    vntr_eval, method, rid = task
    gfa = os.path.join(vntr_eval, "candidates", method, rid + ".gfa")
    if os.path.exists(gfa) and os.path.getsize(gfa) > 0:
        return rid, "cached", 0.0
    sys.path.insert(0, PGREALIGN)
    from pgrealign.realign import MAX_CELLS
    from pgrealign.star import Abpoa, Scores, star_msa
    hap = os.path.join(vntr_eval, "regions", rid, "hap32.fa")
    names, seqs = read_fasta(hap)
    weight = collections.Counter(seqs)
    distinct = list(dict.fromkeys(seqs))
    alleles = [("a%d" % k, s, weight[s]) for k, s in enumerate(distinct)]
    start = time.time()
    try:
        if len(distinct) == 1:
            rows = {"a0": distinct[0]}
        else:
            with Abpoa(exe=ABPOA, scores=Scores.parse("-M 8 -X 16 -O 16,96 -E 8,1")) as align:
                rows = star_msa(alleles, align=align, max_cells=MAX_CELLS).rows
    except Exception as error:
        return rid, "star failed: %s" % error, time.time() - start
    seconds = time.time() - start
    msa = os.path.join(vntr_eval, "work", "plegma", method, rid + ".msa.fa")
    row = {s: rows["a%d" % k] for k, s in enumerate(distinct)}
    with open(msa, "w") as f:
        for n, s in zip(names, seqs):
            f.write(">%s\n%s\n" % (n, row[s]))
    g = subprocess.run(["python3", "-I", os.path.join(vntr_eval, "tools", "msa_graph.py"), msa, hap, gfa, "--quiet"],
                       capture_output=True, text=True)
    return rid, "ok" if g.returncode == 0 else "graph failed: " + g.stderr[-200:], seconds


def init(pgrealign, abpoa):
    global PGREALIGN, ABPOA
    PGREALIGN, ABPOA = pgrealign, abpoa


def main():
    args = sys.argv[1:]
    jobs = 8
    if "--jobs" in args:
        k = args.index("--jobs")
        jobs = int(args[k + 1])
        del args[k:k + 2]
    vntr_eval, pgrealign, abpoa = args[:3]
    method = args[3] if len(args) > 3 else "pgrealign_star"
    os.makedirs(os.path.join(vntr_eval, "candidates", method), exist_ok=True)
    os.makedirs(os.path.join(vntr_eval, "work", "plegma", method), exist_ok=True)
    regions = sorted(d for d in os.listdir(os.path.join(vntr_eval, "regions"))
                     if os.path.exists(os.path.join(vntr_eval, "regions", d, "hap32.fa")))
    with Pool(jobs, initializer=init, initargs=(pgrealign, abpoa)) as pool:
        results = pool.map(one, [(vntr_eval, method, rid) for rid in regions], chunksize=1)
    bad = [(r, m) for r, m, _ in results if m not in ("ok", "cached")]
    print("%d regions, %d failed, %.0f s summed" % (len(results), len(bad), sum(s for _, _, s in results)))
    for r, m in bad:
        print("  %s: %s" % (r, m))


if __name__ == "__main__":
    main()
