#!/usr/bin/env python3
"""vntr-eval's hard set (its 45 chr20 regions whose realignment agrees worst with the pairwise optima) as a fast
loop: plegma's MSAs scored by pgrealign qc's pair F1, the same pairs and optima vntr-eval's strategies are scored
on (tools/hardset_strategies.py's `optima` and `score_region`), without mapping or genotyping.

usage: python3 -I experiments/hardset/hardset.py <vntr-eval> <packages.jsonl.gz> <work dir> <label> [--flags "..."]
           [--regions a,b,... | --pick n] [--jobs n] [--timeout s] [--star] [--gapnorm]

- Regions: `--regions`, or the `--pick` (6) cheapest of the 45 (distinct alleles x the longest's length squared).
- Alleles: each region's, from pgrealign extract's packages (the same graph's); a region's optima are cached in
  <work dir>/optima/<id>.jsonl.gz (qc's 20 random allele pairs and CHM13's allele against the 20 heaviest, aligned
  by pgrealign's star.align_pair).
- Runs `plegma <flags> -r 1` (`-A` by default) on the region's distinct alleles, as experiments/chr20/realign_plegma.py
  does, into <work dir>/msa/<label>/<id>.fa, and with `--star` pgrealign's medoid star too (label `star`), which
  scored on the same alleles and pairs (the hard set's `pair_f1_realigned` comes from vntr-eval's own extract, whose
  alleles differ); with `--gapnorm` each MSA is also scored with its gaps left-normalised (vntr-eval's gap_norm.py,
  its `star_gapnorm`), as F1 counts a co-optimal gap placement as wrong.
- Prints per region the pair F1 of the original graph and of the star (from hardset/chr20.tsv), and this run's pair
  F1, ref F1, cost ratio and seconds; then the medians.

Needs VNTR_ABPOA (abPOA, for the optima and the star) and VNTR_PGREALIGN (pgrealign's checkout) in the environment.
"""
import csv
import gzip
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from multiprocessing import Pool

HS = None


def load(vntr_eval):
    global HS
    if HS is not None:
        return
    sys.path.insert(0, os.path.join(vntr_eval, "tools"))
    import hardset_strategies  # noqa: E402

    HS = hardset_strategies


def regions(vntr_eval):
    with open(os.path.join(vntr_eval, "hardset", "chr20.tsv")) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def packages(path, ids):
    out = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in ids:
                out[p["id"]] = p
                if len(out) == len(ids):
                    break
    return out


def alleles_of(package):
    """[(id, seq, weight)], distinct sequences only (the first allele spelling each), as the hard set's FASTA."""
    seen = set()
    out = []
    for a in package["alleles"]:
        s = a["seq"].upper()
        if s in seen:
            continue
        seen.add(s)
        out.append((a["id"], s, a["weight"]))
    return out


def cache_optima(task):
    vntr_eval, region, alleles, path = task
    load(vntr_eval)
    if not os.path.exists(path):
        with tempfile.TemporaryDirectory() as wd:
            HS.optima(region, alleles, path, wd, 9)
    return region["id"]


def plegma_msa(plegma, flags, alleles, timeout):
    """plegma's rows ({allele id: row}) for the distinct `alleles`, and its seconds."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.fa")
        with open(path, "w") as f:
            for k, (_, s, _) in enumerate(alleles):
                f.write(">s%d\n%s\n" % (k, s))
        start = time.time()
        out = subprocess.run([plegma] + flags + ["-r", "1", path], capture_output=True, text=True, timeout=timeout)
        seconds = time.time() - start
        if out.returncode != 0:
            raise RuntimeError(out.stderr.strip()[-200:])
    rows = {}
    name = None
    for line in out.stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name is not None:
            rows[name].append(line.strip().upper())
    return {a[0]: "".join(rows["s%d" % k]) for k, a in enumerate(alleles)}, seconds


def star_msa(region, alleles):
    with tempfile.TemporaryDirectory() as wd:
        start = time.time()
        rows = HS.s_star(alleles, wd, region["ref_allele"])
        return rows, time.time() - start


def one(task):
    vntr_eval, label, plegma, flags, region, alleles, work, timeout, gapnorm = task
    load(vntr_eval)
    rid = region["id"]
    out = os.path.join(work, "msa", label, rid + ".fa")
    try:
        if label == "star":
            rows, seconds = star_msa(region, alleles)
            if isinstance(rows, list):
                rows = {a[0]: r for a, r in zip(alleles, rows)}
        else:
            rows, seconds = plegma_msa(plegma, flags, alleles, timeout)
    except subprocess.TimeoutExpired:
        return rid, {"status": "timeout"}
    except Exception as error:  # noqa: BLE001
        return rid, {"status": "failed: %s" % str(error)[-120:]}
    for a, s, _ in alleles:
        if rows[a].replace("-", "") != s:
            return rid, {"status": "row %s does not spell its allele" % a}
    with open(out, "w") as f:
        for a, _, _ in alleles:
            f.write(">%s\n%s\n" % (a, rows[a]))
    with gzip.open(os.path.join(work, "optima", rid + ".jsonl.gz"), "rt") as f:
        opt = [json.loads(line) for line in f]
    out = {"status": "ok", "seconds": seconds}
    out.update(scores(rows, opt))
    if gapnorm:
        # The same MSA, every row's gap runs shifted left against the column consensus (vntr-eval's gap_norm.py).
        import gap_norm  # noqa: E402

        names = [a for a, _, _ in alleles]
        moved, _ = gap_norm.normalise([rows[a] for a in names], "left")
        out["gapnorm"] = scores(dict(zip(names, moved)), opt)
    return rid, out


def scores(rows, opt):
    sc = HS.score_region(rows, opt)
    return {
        "pair_f1": HS.f1(sc["pair_tp"], sc["pair_opt"], sc["pair_msa"]),
        "ref_f1": HS.f1(sc["ref_tp"], sc["ref_opt"], sc["ref_msa"]) if sc["ref_opt"] else None,
        "cost_ratio": sc["pair_cost"] / sc["pair_cost_opt"] if sc["pair_cost_opt"] else None,
    }


def main():
    args = sys.argv[1:]

    def option(name, default):
        if name in args:
            k = args.index(name)
            value = args[k + 1]
            del args[k : k + 2]
            return value
        return default

    flags = option("--flags", "-A").split()
    pick = int(option("--pick", "6"))
    chosen = option("--regions", "")
    jobs = int(option("--jobs", "6"))
    timeout = int(option("--timeout", "1800"))
    with_star = "--star" in args
    gapnorm = "--gapnorm" in args
    args = [a for a in args if a not in ("--star", "--gapnorm")]
    vntr_eval, package_path, work, label = args[:4]
    plegma = os.environ.get("PLEGMA", os.path.join(os.path.dirname(__file__), "..", "..", "build", "plegma"))
    load(vntr_eval)
    every = regions(vntr_eval)
    if chosen:
        wanted = chosen.split(",")
        picked = [r for r in every if r["id"] in wanted]
    else:
        picked = sorted(every, key=lambda r: int(r["alleles"]) * int(r["longest"]) ** 2)[:pick]
    found = packages(package_path, {r["id"] for r in picked})
    picked = [r for r in picked if r["id"] in found]
    alleles = {r["id"]: alleles_of(found[r["id"]]) for r in picked}
    os.makedirs(os.path.join(work, "optima"), exist_ok=True)
    labels = [label] + (["star"] if with_star else [])
    for lab in labels:
        os.makedirs(os.path.join(work, "msa", lab), exist_ok=True)
    with Pool(jobs) as pool:
        tasks = [(vntr_eval, r, alleles[r["id"]], os.path.join(work, "optima", r["id"] + ".jsonl.gz")) for r in picked]
        list(pool.imap_unordered(cache_optima, tasks))
        results = {}
        tasks = [(vntr_eval, lab, plegma, flags, r, alleles[r["id"]], work, timeout, gapnorm) for lab in labels
                 for r in picked]
        for (rid, res), task in zip(pool.imap(one, tasks), tasks):
            results[(task[1], rid)] = res
    shown = [(lab, False) for lab in labels] + ([(lab, True) for lab in labels] if gapnorm else [])
    names = [lab + ("+gapnorm" if g else "") for lab, g in shown]
    print("%-20s %-9s %7s %6s | %6s %6s | %s" % ("region", "class", "alleles", "graph", "star*", "", " ".join(
        "%-26s" % name for name in names)))
    medians = {name: [] for name in names}
    for r in picked:
        cells = []
        for (lab, g), name in zip(shown, names):
            res = results[(lab, r["id"])]
            if res["status"] == "ok":
                sc = res["gapnorm"] if g else res
                medians[name].append(sc["pair_f1"])
                cells.append("F1 %.4f cost %.3f %6.1fs" % (sc["pair_f1"], sc["cost_ratio"] or 0, res["seconds"]))
            else:
                cells.append("%-26s" % res["status"][:26])
        print("%-20s %-9s %7s %6s | %6s %6s | %s" % (r["id"], r["classes"][:9], len(alleles[r["id"]]),
                                                   r["pair_f1_graph"], r["pair_f1_realigned"], "", " ".join(cells)))
    print("median pair F1: graph %.4f, star* %.4f | %s" % (
        statistics.median(float(r["pair_f1_graph"]) for r in picked),
        statistics.median(float(r["pair_f1_realigned"]) for r in picked),
        ", ".join("%s %.4f (%d)" % (lab, statistics.median(v), len(v)) for lab, v in medians.items() if v)))
    print("(star*: the hard set's pair_f1_realigned, pgrealign's medoid star on vntr-eval's own extract)")


if __name__ == "__main__":
    main()
