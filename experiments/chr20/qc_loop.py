#!/usr/bin/env python3
"""A fast loop on a fixed subset of pgrealign's chr20 regions: each MSA set scored by pgrealign qc, no graphs, no
genotyping.

usage: python3 -I experiments/chr20/qc_loop.py <work dir> <regions file> <out dir> <label>=<msas.jsonl.gz | plegma
           flags> ... [--scores asm5] [--jobs n] [--graph graph.gbz --vg vg] [--deployed label] [--coopt]

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

With --graph (the chr20 graph extract read), qc measures the graph's own alignment of each region's pairs too (its
graph_* columns, for every label's regions), and with `--deployed label` that label is also shown as pgrealign
deploys it: its MSA where realign made one, the graph's own alignment where it refused the region (not_homologous,
timeout, ...), over every region of the subset.

qc counts a letter pair right only when the pair's one optimal alignment it computes (left-normalised) aligns it, so
an MSA that places an indel at another place the same cost allows scores lower. With --coopt each label's pairs
are scored again counting a letter pair right when some optimal alignment aligns it (build/anchors --coopt: forward
and backward DP under the scores' two-piece costs; ANCHORS): the co-optimal pair and reference F1, which no choice
among equally good placements changes.

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
# The same scores as two-piece costs (a match 0): mismatch 2 (M + X), gap opens 2 O, extensions 2 E + M.
TWO_PIECE = {"default": "48,32,24,192,10", "asm5": "40,78,7,162,3"}
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


def coopt(msas, packages, scores, jobs, workdir):
    """Co-optimal counts of a label's MSAs: {"pair": [on], "ref": [on]} summed over its regions' pairs, `on` the MSA's
    letter pairs some optimal alignment aligns. Each pair's co-optimal letter pairs (build/anchors --coopt-set) are
    kept on disk in <out dir>/coopt, keyed by the scores and the two sequences: they do not depend on the MSA."""
    import tempfile

    rows_of = {}
    with gzip.open(msas, "rt") as f:
        for line in f:
            m = json.loads(line)
            if m.get("status") == "ok" and len(m["rows"]) > 1:
                rows_of[m["id"]] = m["rows"]
    cache = os.path.join(workdir, "coopt")
    todo = {}
    work = []
    with gzip.open(packages, "rt") as f:
        for line in f:
            p = json.loads(line)
            rows = rows_of.get(p["id"])
            if rows is None or len(p["alleles"]) < 2:
                continue
            seqs = {a["id"]: a["seq"].upper() for a in p["alleles"]}
            for kind, a, b in qc.compared_pairs(p, rows, 20, 20):
                key = hashlib.sha1(("%s|%s|%s" % (TWO_PIECE[scores], seqs[a], seqs[b])).encode()).hexdigest()
                path = os.path.join(cache, key[:2], key + ".txt")
                if not os.path.exists(path):
                    todo[key] = (seqs[a], seqs[b], path)
                work.append((kind, path, qc.induced(rows[a], rows[b])))
    if todo:
        anchors = os.environ.get("ANCHORS", os.path.expanduser("~/plegma-next/build/anchors"))
        items = sorted(todo.values(), key=lambda t: -len(t[0]) * len(t[1]))
        chunks = [items[k::jobs] for k in range(jobs)]
        with tempfile.TemporaryDirectory(dir=workdir) as tmp:
            procs = []
            for k, chunk in enumerate(chunks):
                if not chunk:
                    continue
                path = os.path.join(tmp, "%d.tsv" % k)
                with open(path, "w") as f:
                    f.writelines("%s\t%s\n" % (a, b) for a, b, _ in chunk)
                procs.append((chunk, subprocess.Popen([anchors, "--coopt-set", path, "--two-piece",
                                                       TWO_PIECE[scores]], stdout=subprocess.PIPE, text=True)))
            for chunk, proc in procs:
                for (_, _, out_path), row in zip(chunk, proc.stdout):
                    os.makedirs(os.path.dirname(out_path), exist_ok=True)
                    with open(out_path + ".tmp", "w") as f:
                        f.write(row)
                    os.replace(out_path + ".tmp", out_path)
                proc.wait()
    totals = {"pair": [0], "ref": [0]}
    parsed = {}
    for kind, path, (ra, rb) in work:
        if path not in parsed:
            parsed[path] = parse_coopt(open(path).read())
        totals[kind][0] += on_optimum(parsed[path], ra, rb)
    return totals


def parse_coopt(line):
    """A --coopt-set line as (lo, hi, extra): each letter of a's first interval of co-optimal partners (lo > hi for
    none) and its further intervals, {letter: [(lo, hi)]}, for the few with more than one."""
    import numpy as np

    tokens = line.split(" ", 2)[2].rstrip("\n").split(",")
    lo = np.empty(len(tokens), dtype=np.int64)
    hi = np.empty(len(tokens), dtype=np.int64)
    extra = {}
    for i, t in enumerate(tokens):
        if t == ".":
            lo[i], hi[i] = 1, 0
            continue
        intervals = t.split(":")
        first = intervals[0].partition("-")
        lo[i], hi[i] = int(first[0]), int(first[2] or first[0])
        if len(intervals) > 1:
            extra[i] = [(int(x.partition("-")[0]), int(x.partition("-")[2] or x.partition("-")[0])) for x in intervals[1:]]
    return lo, hi, extra


def on_optimum(sets, ra, rb):
    """How many letter pairs the rows `ra`, `rb` align that some optimal alignment aligns."""
    import numpy as np

    lo, hi, extra = sets
    a = np.frombuffer(ra.encode(), dtype=np.uint8) != ord("-")
    b = np.frombuffer(rb.encode(), dtype=np.uint8) != ord("-")
    both = a & b
    i = (np.cumsum(a) - 1)[both]
    j = (np.cumsum(b) - 1)[both]
    hit = (j >= lo[i]) & (j <= hi[i])
    on = int(hit.sum())
    if extra:
        for k in np.nonzero(~hit)[0]:
            for x, y in extra.get(int(i[k]), ()):
                if x <= j[k] <= y:
                    on += 1
                    break
    return on


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
    graph = option(args, "--graph", "")
    vg = option(args, "--vg", os.environ.get("VG", "vg"))
    deployed = option(args, "--deployed", "")
    with_coopt = "--coopt" in args
    args = [a for a in args if a != "--coopt"]
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
                       scores=SCORES[s], workdir=out, graph=graph or None, vg=vg)
            if with_coopt:
                on = coopt(msas, sub, s, jobs, out)
                for kind in ("pair", "ref"):
                    total = d[kind + "_opt"] + d[kind + "_msa"]
                    d["co_" + kind + "_f1"] = round(2 * on[kind][0] / total, 4) if total else None
            rows.append((label, s, d, made, time.time() - t))
    print("%d regions" % len(ids))
    header = "%-10s %-8s %8s %8s %13s %13s %7s %8s" % ("label", "scores", "pair F1", "ref F1", "ref_sv/opt",
                                                        "forced/opt", "growth", "regions")
    print(header + " %9s" % "seconds")
    for label, s, d, made, took in rows:
        print(line(label, s, d) + " %4.0f+%-4.0f" % (made, took))
    if with_coopt:
        print("co-optimal (a letter pair right when some optimal alignment aligns it):")
        for label, s, d, _, _ in rows:
            print("%-10s %-8s co pair F1 %.4f  co ref F1 %.4f" % (label, s, d["co_pair_f1"], d["co_ref_f1"]))
    # The same over the regions every label measured.
    per = {}
    for label, s, _, _, _ in rows:
        with open(os.path.join(out, "%s.%s.tsv" % (label, s))) as f:
            names = f.readline().rstrip("\n").split("\t")
            per[label, s] = {r[0]: dict(zip(names, [r[0]] + [int(x) for x in r[1:]]))
                             for r in (x.rstrip("\n").split("\t") for x in f)}
    if graph:
        print("the graph's own alignment (qc --graph), over each label's regions:")
        for (label, s), regions in per.items():
            d = qc.summarise(list(regions.values()))
            print("%-10s %-8s graph pair F1 %.4f, graph ref F1 %.4f, %d regions" % (
                label, s, d["graph_pair_f1"], d["graph_ref_f1"], d["regions"]))
    if deployed and graph:
        # The label as pgrealign deploys it: its MSA's counts where it made one, the graph's own where it refused.
        for s in scores:
            mine = per.get((deployed, s), {})
            full = max((v for (lab, sc), v in per.items() if sc == s), key=len)
            rows = []
            for rid, r in full.items():
                if rid in mine:
                    rows.append(mine[rid])
                else:
                    g = dict(r)
                    for k in ("pair_tp", "pair_msa", "ref_tp", "ref_msa"):
                        g[k] = r["graph_" + k]
                    g["sequence"] = r["graph_sequence"]
                    for k in ("ref_sv", "forced"):
                        g[k] = 0
                    rows.append(g)
            d = qc.summarise(rows)
            print("%s as deployed (%d MSAs, %d regions at the graph's own alignment), %s: pair F1 %.4f, ref F1 %.4f"
                  % (deployed, len(mine), len(full) - len(mine), s, d["pair_f1"], d["ref_f1"]))
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
