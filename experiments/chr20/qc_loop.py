#!/usr/bin/env python3
"""A fast loop on a fixed subset of pgrealign's chr20 regions: each MSA set scored by pgrealign qc, no graphs, no
genotyping.

usage: python3 -I experiments/chr20/qc_loop.py <work dir> <regions file> <out dir> <label>=<msas.jsonl.gz | plegma
           flags> ... [--scores asm5] [--jobs n]

For each label, its MSAs over the regions in <regions file> (one id per line): a msas.jsonl.gz (the star's, an
arm's) cut to them, or plegma's for those flags (realign_plegma.py on those regions' packages, the CHM13 flanks
given for -F; PLEGMA and KARAN the binaries). pgrealign qc then measures them under each `--scores`: `asm5`
(minimap2's, -A's costs: -M 1 -X 19 -O 39,81 -E 3,1, the default, as qc asks for the scores an MSA is built with) or
`default` (pgrealign's, abPOA's -M 8 -X 16 -O 16,96 -E 8,1, the star's):
qc compares the alignment the MSA implies for 20 random pairs of alleles, and for CHM13's allele against the 20
heaviest, with each pair's optimal alignment under those scores.

qc runs here in-process with its pairs' optimal alignments (pgrealign.qc.align_pair) cached on disk in
<out dir>/optima, keyed by the scores and the two sequences: they do not depend on the MSA, so only a subset's first
run pays for them, and every later label and iteration scores in seconds.

Prints per label and scores the pair and reference F1, the long gap runs (50 or more gap columns) against the
optimal alignments', the forced-homology columns against theirs, the growth (sequence over the longest allele) and
the regions measured, then the same over the regions every label measured (the star lacks a few); writes
<out dir>/<label>.<scores>.tsv per region. Needs pgrealign (PGREALIGN,
~/rivals-src/pgrealign) and abPOA (ABPOA, ~/plegma-bench/msa/bin/abpoa).
"""
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("PGREALIGN", os.path.expanduser("~/rivals-src/pgrealign")))
from pgrealign import qc  # noqa: E402
from pgrealign.star import Scores  # noqa: E402

SCORES = {"default": Scores(match=8, mismatch=16, gap_open=(16, 96), gap_ext=(8, 1)),
          "asm5": Scores(match=1, mismatch=19, gap_open=(39, 81), gap_ext=(3, 1))}
CACHE = os.environ.get("QC_OPTIMA_CACHE", "")
_align_pair = qc.align_pair


def cached_align_pair(centre, seq, align, banded=False):
    """pgrealign's align_pair, kept on disk by (scores, centre, seq, banded)."""
    if not CACHE:
        return _align_pair(centre, seq, align, banded=banded)
    flags = " ".join(align.scores.flags())
    key = hashlib.sha1(("%s|%d|%s|%s" % (flags, banded, centre, seq)).encode()).hexdigest()
    path = os.path.join(CACHE, key[:2], key + ".json")
    if os.path.exists(path):
        return tuple(json.load(open(path)))
    rows = _align_pair(centre, seq, align, banded=banded)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".%d" % os.getpid(), "w") as f:
        json.dump(rows, f)
    os.replace(path + ".%d" % os.getpid(), path)
    return rows


# At import, so qc's worker processes (spawned, forkserver or forked) align through the cache too.
qc.align_pair = cached_align_pair


def option(args, name, default):
    if name in args:
        i = args.index(name)
        value = args[i + 1]
        del args[i : i + 2]
        return value
    return default


def subset(path, ids, out):
    with gzip.open(path, "rt") as f, gzip.open(out, "wt") as o:
        for line in f:
            if json.loads(line)["id"] in ids:
                o.write(line)


def main():
    global CACHE
    args = sys.argv[1:]
    scores = option(args, "--scores", "asm5").split(",")
    jobs = int(option(args, "--jobs", "8"))
    work, regions, out = args[:3]
    labels = [x.split("=", 1) for x in args[3:]]
    os.makedirs(out, exist_ok=True)
    CACHE = os.environ["QC_OPTIMA_CACHE"] = os.path.abspath(os.path.join(out, "optima"))
    ids = set(open(regions).read().split())
    extract = os.path.join(work, "star", "chr20", "extract")
    # The packages of the regions, kept per regions file (a smaller subset shares the directory and its optima).
    sub = os.path.join(out, "packages.%s.jsonl.gz" % os.path.splitext(os.path.basename(regions))[0])
    if not os.path.exists(sub):
        subset(os.path.join(extract, "packages.jsonl.gz"), ids, sub)
    abpoa = os.environ.get("ABPOA", os.path.expanduser("~/plegma-bench/msa/bin/abpoa"))
    rows = []
    for label, source in labels:
        started = time.time()
        msas = os.path.join(out, label + ".msas.jsonl.gz")
        if source.endswith(".jsonl.gz"):
            subset(source, ids, msas)
        else:
            if os.path.exists(msas):
                os.remove(msas)
            subprocess.run(["python3", "-I", os.path.join(HERE, "realign_plegma.py"), os.environ["PLEGMA"],
                            os.environ["KARAN"], sub, msas, "--jobs", str(jobs), "--timeout", "600", "--flags",
                            source, "--reference", os.path.join(work, "calls", "chr20.fa")], check=True,
                           stdout=subprocess.DEVNULL, stderr=open(os.path.join(out, label + ".realign.err"), "w"))
        made = time.time() - started
        for s in scores:
            t = time.time()
            d = qc.run(extract, msas, os.path.join(out, "%s.%s.tsv" % (label, s)), jobs=jobs, abpoa=abpoa,
                       scores=SCORES[s], workdir=out)
            rows.append((label, s, d, made, time.time() - t))
    print("%d regions" % len(ids))
    header = "%-10s %-8s %8s %8s %13s %13s %7s %8s" % ("label", "scores", "pair F1", "ref F1", "ref_sv/opt",
                                                        "forced/opt", "growth", "regions")
    print(header + " %9s" % "seconds")
    for label, s, d, made, took in rows:
        print(line(label, s, d) + " %4.0f+%-4.0f" % (made, took))
    # The same over the regions every label measured.
    per = {}
    for label, s, _, _, _ in rows:
        with open(os.path.join(out, "%s.%s.tsv" % (label, s))) as f:
            names = f.readline().rstrip("\n").split("\t")
            per[label, s] = {r[0]: dict(zip(names, [r[0]] + [int(x) for x in r[1:]]))
                             for r in (x.rstrip("\n").split("\t") for x in f)}
    common = set.intersection(*(set(v) for v in per.values())) if per else set()
    if any(len(v) != len(common) for v in per.values()):
        print("over the %d regions every label measured:" % len(common))
        print(header)
        for (label, s), regions in per.items():
            print(line(label, s, qc.summarise([regions[r] for r in sorted(common)])))


def line(label, s, d):
    return "%-10s %-8s %8.4f %8.4f %13s %13s %7.4f %8d" % (
        label, s, d["pair_f1"], d["ref_f1"], "%d/%d" % (d["ref_sv"], d["ref_sv_opt"]),
        "%d/%d" % (d["forced"], d["forced_opt"]), d["growth"], d["regions"])


if __name__ == "__main__":
    main()
