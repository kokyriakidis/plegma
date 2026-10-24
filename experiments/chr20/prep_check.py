#!/usr/bin/env python3
"""Checks of a local_call.py prep, region by region: what every later step takes from it.

usage: python3 -I experiments/chr20/prep_check.py <prep dir> <packages.jsonl.gz> [--regions file] [--jobs n]

Per region: the window (chr20.fa) spelling left flank + CHM13's allele + right flank; R1 and R2 paired name for name;
the truth's small and structural records inside the interior (the part the MSA aligns, which is what is scored), in
the flanks, and reaching past the window; the share of the interior the benchmark beds cover; the reads, and their
depth over the left flank, the interior and the right flank (median count of K-mers the window holds once). Writes
<prep dir>/prep_qc.tsv and prints the regions failing a check: the window not spelling, mates out of step, a flank
under MIN_FLANK_DEPTH, an interior under STARVED of the flanks' depth (reads the input lacks: reads.sh takes GRCh38's
chr20 alone), an interior the benchmark covers under MIN_BED.
"""
import gzip
import json
import os
import statistics
import sys
from multiprocessing import Pool

K = 31
COMP = str.maketrans("ACGT", "TGCA")
MIN_FLANK_DEPTH = 10
STARVED = 0.4
MIN_BED = 0.5


def canonical(seq):
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    for i in range(n - K + 1):
        a = seq[i : i + K]
        b = rc[n - K - i : n - i]
        yield a if a <= b else b


def names(path):
    with open(path) as f:
        return [line[1:].split()[0].split("/")[0] for k, line in enumerate(f) if k % 4 == 0]


def truth_counts(path, lo, hi, length):
    inside = flank = past = 0
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            w = line.split("\t", 4)
            start, end = int(w[1]) - 1, int(w[1]) - 1 + len(w[3])
            if end > length:
                past += 1
            elif lo <= start and end <= hi:
                inside += 1
            else:
                flank += 1
    return inside, flank, past


def check(task):
    prep_dir, pkg = task
    rid = pkg["id"]
    d = os.path.join(prep_dir, rid)
    meta = json.load(open(os.path.join(d, "region.json")))
    window = "".join(x.strip() for x in open(os.path.join(d, "chr20.fa")) if not x.startswith(">")).upper()
    ref = next((a["seq"].upper() for a in pkg["alleles"] if a["id"] == pkg.get("ref")), None)
    lo, hi = len(meta["left"]), len(window) - len(meta["right"])
    out = {"region": rid, "window": len(window), "interior": hi - lo,
           "window_spells": int(ref is not None and meta["left"] + ref + meta["right"] == window)}
    r1, r2 = names(os.path.join(d, "R1.fq")), names(os.path.join(d, "R2.fq"))
    out["pairs"], out["singles"] = len(r1), len(names(os.path.join(d, "single.fq")))
    out["mates_in_step"] = int(r1 == r2)
    for kind in ("smvar", "stvar"):
        i, f, p = truth_counts(os.path.join(d, "truth.chr20.%s.vcf.gz" % kind), lo, hi, len(window))
        out[kind + "_interior"], out[kind + "_flanks"], out[kind + "_past"] = i, f, p
    covered = 0
    for line in open(os.path.join(d, "truth.chr20.stvar.bed")):
        w = line.split("\t")
        covered += max(0, min(int(w[2]), hi) - max(int(w[1]), lo))
    out["bed_interior"] = round(covered / max(1, hi - lo), 3)
    parts = {"left": window[:lo], "interior": window[lo:hi], "right": window[hi:]}
    once = {}
    counts = {}
    for km in canonical(window):
        counts[km] = counts.get(km, 0) + 1
    for p, s in parts.items():
        once[p] = [km for km in set(canonical(s)) if counts.get(km) == 1]
    seen = dict.fromkeys((km for v in once.values() for km in v), 0)
    for name in ("R1", "R2", "single"):
        with open(os.path.join(d, name + ".fq")) as f:
            for k, line in enumerate(f):
                if k % 4 == 1:
                    for km in canonical(line.strip().upper()):
                        if km in seen:
                            seen[km] += 1
    for p, v in once.items():
        out["depth_" + p] = statistics.median(seen[km] for km in v) if v else ""
    flags = []
    if not out["window_spells"]:
        flags.append("window")
    if not out["mates_in_step"]:
        flags.append("mates")
    flank = [x for x in (out["depth_left"], out["depth_right"]) if x != ""]
    if flank and min(flank) < MIN_FLANK_DEPTH:
        flags.append("flank_depth")
    if flank and out["depth_interior"] != "" and out["depth_interior"] < STARVED * statistics.mean(flank):
        flags.append("starved")
    if out["bed_interior"] < MIN_BED:
        flags.append("bed")
    out["flags"] = ",".join(flags)
    return out


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    regions = None
    jobs = 8
    if "--regions" in args:
        regions = set(open(args[args.index("--regions") + 1]).read().split())
        del args[args.index("--regions") : args.index("--regions") + 2]
    if "--jobs" in args:
        jobs = int(args[args.index("--jobs") + 1])
        del args[args.index("--jobs") : args.index("--jobs") + 2]
    prep_dir, packages = args
    pkgs = []
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if os.path.exists(os.path.join(prep_dir, p["id"], "region.json")) and (regions is None or p["id"] in regions):
                pkgs.append(p)
    with Pool(jobs) as pool:
        rows = list(pool.imap_unordered(check, [(prep_dir, p) for p in pkgs]))
    rows.sort(key=lambda r: r["region"])
    cols = list(rows[0])
    with open(os.path.join(prep_dir, "prep_qc.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("%d regions" % len(rows))
    for flag in ("window", "mates", "flank_depth", "starved", "bed"):
        hit = [r["region"] for r in rows if flag in r["flags"].split(",")]
        print("  %-12s %3d  %s" % (flag, len(hit), " ".join(hit[:8])))
    print("  truth records in the flanks (not scored): small %d, structural %d; past the window: %d, %d" % (
        sum(r["smvar_flanks"] for r in rows), sum(r["stvar_flanks"] for r in rows),
        sum(r["smvar_past"] for r in rows), sum(r["stvar_past"] for r in rows)))
    for p in ("left", "interior", "right"):
        v = [r["depth_" + p] for r in rows if r["depth_" + p] != ""]
        print("  depth %-9s median %.0f, 10%% %.0f" % (p, statistics.median(v), sorted(v)[len(v) // 10]))


if __name__ == "__main__":
    main()
