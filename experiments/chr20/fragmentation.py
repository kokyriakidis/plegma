#!/usr/bin/env python3
"""Whether an MSA keeps a large insertion or deletion against CHM13 as one event (vg call reports variants against
CHM13's path, so an event split into pieces is a missed truth variant and spurious calls): a fast check, region by
region, without mapping or calling.

usage: python3 -I experiments/chr20/fragmentation.py <packages.jsonl.gz> <regions a,b,...|file> <out dir>
           <label>=<msas.jsonl.gz | plegma flags> [...] [--jobs n]

For each region, CHM13's allele against every other allele (each weighing its haplotypes): the gaps of MIN letters or
more in their alignment, the optimal pairwise one (pgrealign's star.align_pair under its scores, the alignment qc
compares against) and each label's induced one (a label is a pgrealign msas.jsonl.gz, or `plegma` flags run on the
region's distinct alleles). An optimal gap is `kept` when one induced gap on its side covers at least KEPT of it, else
`split`; induced gaps that no optimal gap overlaps are `extra`. Prints per region and label kept / split / extra
(weighted) and the totals. Needs VNTR_ABPOA, VNTR_PGREALIGN and, for flags, PLEGMA.
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile
from multiprocessing import Pool

MIN = 50
KEPT = 0.8


def gaps(ra, rb):
    """The gap runs of MIN letters or more in a two-row alignment: (side, start in CHM13's letters, length); side `D`
    letters of the first row (CHM13's) unaligned (a deletion in the other), `I` letters of the second."""
    out = []
    i = 0
    run_side, run_start, run_len = None, 0, 0
    for x, y in zip(ra, rb):
        if x == "-" and y == "-":
            continue
        side = "D" if y == "-" else "I" if x == "-" else None
        if side != run_side:
            if run_side and run_len >= MIN:
                out.append((run_side, run_start, run_len))
            run_side, run_start, run_len = side, i, 0
        if side:
            run_len += 1
        if x != "-":
            i += 1
    if run_side and run_len >= MIN:
        out.append((run_side, run_start, run_len))
    return out


def compare(opt, mine):
    kept = split = 0
    for side, start, length in opt:
        if any(s == side and l >= KEPT * length and abs(st - start) <= max(length, 50) for s, st, l in mine):
            kept += 1
        else:
            split += 1
    extra = sum(1 for s, st, l in mine if not any(o[0] == s and st < o[1] + o[2] + 50 and o[1] < st + l + 50
                                                  for o in opt))
    return kept, split, extra


def induced(ra, rb):
    keep = [i for i in range(len(ra)) if ra[i] != "-" or rb[i] != "-"]
    return "".join(ra[i] for i in keep), "".join(rb[i] for i in keep)


def plegma_rows(flags, seqs, copies=None):
    """plegma's rows for `seqs` (CHM13's first); CHM13_WEIGHT=k feeds CHM13's allele k times over, so -A's objective,
    which weighs each distinct sequence by its copies, counts CHM13's pairs k times."""
    extra = int(os.environ.get("CHM13_WEIGHT", "1"))
    weights = copies
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "in.fa")
        with open(path, "w") as f:
            for k, s in enumerate(seqs):
                f.write(">s%d\n%s\n" % (k, s))
                # WEIGHTED=1: once per haplotype spelling it, so -A weighs each sequence as the panel does.
                for c in range((weights or {}).get(s, 1) - 1):
                    f.write(">s%d_%d\n%s\n" % (k, c, s))
            for k in range(extra - 1):
                f.write(">extra%d\n%s\n" % (k, seqs[0]))
        out = subprocess.run([os.environ["PLEGMA"]] + flags.split() + ["-r", "1", path], capture_output=True,
                             text=True, check=True).stdout
    rows, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            rows[name] = []
        elif name is not None:
            rows[name].append(line.strip().upper())
    return ["".join(rows["s%d" % k]) for k in range(len(seqs))]


def optimal_gaps(pkg, ref, cache):
    """{allele id: the optimal alignment's gaps against CHM13's allele}, from `cache` (a JSON file) or computed."""
    if os.path.exists(cache):
        return json.load(open(cache))
    sys.path.insert(0, os.environ["VNTR_PGREALIGN"])
    from pgrealign import star  # noqa: E402

    out = {}
    with tempfile.TemporaryDirectory() as wd, star.Abpoa(os.environ["VNTR_ABPOA"], star.SCORES, workdir=wd) as align:
        for a in pkg["alleles"]:
            if a["id"] == ref["id"] or a["seq"].upper() == ref["seq"].upper():
                continue
            oa, ob = star.align_pair(ref["seq"].upper(), a["seq"].upper(), align,
                                     banded=len(ref["seq"]) * len(a["seq"]) > 4e8)
            out[a["id"]] = gaps(oa, ob)
    with open(cache + ".tmp", "w") as f:
        json.dump(out, f)
    os.replace(cache + ".tmp", cache)
    return out


def cached(task):
    pkg, cache_dir = task
    ref = next((a for a in pkg["alleles"] if a["id"] == pkg.get("ref")), None)
    if ref is not None:
        optimal_gaps(pkg, ref, os.path.join(cache_dir, pkg["id"] + ".json"))
    return pkg["id"]


def one(task):
    pkg, labels, msas, cache_dir = task

    alleles = pkg["alleles"]
    ref = next((a for a in alleles if a["id"] == pkg.get("ref")), None)
    if ref is None:
        return pkg["id"], None
    # CHM13's allele first (PLEGMA_CENTRE=0 centres a star on it).
    distinct = list(dict.fromkeys([ref["seq"].upper()] + [a["seq"].upper() for a in alleles]))
    rows_of = {}
    for label, source in labels:
        if source.endswith(".jsonl.gz"):
            rec = msas.get((label, pkg["id"]))
            if rec is None or rec.get("status") != "ok":
                continue
            rows_of[label] = rec["rows"]
        else:
            weights = None
            if os.environ.get("WEIGHTED") == "1":
                weights = {}
                for a in alleles:
                    weights[a["seq"].upper()] = weights.get(a["seq"].upper(), 0) + a["weight"]
            try:
                rows = plegma_rows(source, distinct, weights)
            except subprocess.CalledProcessError:
                continue
            row_of = dict(zip(distinct, rows))
            rows_of[label] = {a["id"]: row_of[a["seq"].upper()] for a in alleles}
    out = {label: [0, 0, 0] for label in rows_of}
    events = 0
    optimal = optimal_gaps(pkg, ref, os.path.join(cache_dir, pkg["id"] + ".json"))
    for a in alleles:
        if a["id"] not in optimal:
            continue
        opt = [tuple(g) for g in optimal[a["id"]]]
        events += a["weight"] * len(opt)
        for label, rows in rows_of.items():
            if a["id"] not in rows or ref["id"] not in rows:
                continue
            kept, split, extra = compare(opt, gaps(*induced(rows[ref["id"]], rows[a["id"]])))
            out[label][0] += a["weight"] * kept
            out[label][1] += a["weight"] * split
            out[label][2] += a["weight"] * extra
    return pkg["id"], (events, out)


def main():
    args = sys.argv[1:]
    jobs = 8
    if "--jobs" in args:
        jobs = int(args[args.index("--jobs") + 1])
        del args[args.index("--jobs") : args.index("--jobs") + 2]
    packages, chosen, out_dir = args[:3]
    labels = [tuple(x.split("=", 1)) for x in args[3:]]
    ids = open(chosen).read().split() if os.path.exists(chosen) else chosen.split(",")
    wanted = set(ids)
    pkgs = {}
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            if p["id"] in wanted:
                pkgs[p["id"]] = p
    msas = {}
    for label, source in labels:
        if source.endswith(".jsonl.gz"):
            with gzip.open(source, "rt") as f:
                for line in f:
                    r = json.loads(line)
                    if r["id"] in wanted:
                        msas[(label, r["id"])] = r
    os.makedirs(out_dir, exist_ok=True)
    cache_dir = os.environ.get("FRAG_CACHE", os.path.join(out_dir, "optima"))
    os.makedirs(cache_dir, exist_ok=True)
    # One (region, label) per task, so the slowest single run bounds the wall time, not a region's labels in turn;
    # each region's optima first, once.
    with Pool(jobs) as pool:
        list(pool.imap_unordered(cached, [(pkgs[i], cache_dir) for i in ids if i in pkgs]))
    tasks = [(pkgs[i], [lab], {k: v for k, v in msas.items() if k[1] == i}, cache_dir) for i in ids if i in pkgs
             for lab in labels]
    totals = {label: [0, 0, 0] for label, _ in labels}
    merged = {}
    with Pool(jobs) as pool:
        for rid, res in pool.imap_unordered(one, tasks):
            if res is None:
                continue
            events, out = res
            got = merged.setdefault(rid, [events, {}])
            got[1].update(out)
    with open(os.path.join(out_dir, "fragmentation.tsv"), "w") as f:
        f.write("region\tevents\t" + "\t".join("%s_kept\t%s_split\t%s_extra" % (l, l, l) for l, _ in labels) + "\n")
        print("%-22s %7s | %s" % ("region", "events", " | ".join("%-24s" % l for l, _ in labels)))
        for rid in ids:
            if rid not in merged:
                continue
            events, out = merged[rid]
            cells = []
            for label, _ in labels:
                k, s, e = out.get(label, [0, 0, 0])
                for t in range(3):
                    totals[label][t] += (k, s, e)[t]
                cells.append("kept %4d split %4d extra %4d" % (k, s, e))
            f.write("%s\t%d\t%s\n" % (rid, events, "\t".join("%d\t%d\t%d" % tuple(out.get(l, [0, 0, 0]))
                                                              for l, _ in labels)))
            print("%-22s %7d | %s" % (rid, events, " | ".join(cells)), flush=True)
    for label, (k, s, e) in totals.items():
        print("%-12s large events against CHM13 kept whole %d, split %d (%.1f%% kept), extra gaps %d" % (
            label, k, s, 100.0 * k / max(k + s, 1), e))


if __name__ == "__main__":
    main()
