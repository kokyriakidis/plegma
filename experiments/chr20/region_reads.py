#!/usr/bin/env python3
"""The reads each chr20 arm aligns into the tandem-repeat regions, and how well, region by region.

usage: python3 -I experiments/chr20/region_reads.py <work dir> <packages.jsonl.gz> <vg> <arm> <arm> [<arm> ...]

A region's nodes, in every arm's graph: its anchors L and R (pgrealign extract's) keep their sequence and their
place on CHM13's path, and replace numbers densely (the nodes that stay keep their order, each region's new nodes
right after its left anchor), so the region is the node ids from L to R. The first arm is the graph the extract was
made from: there L and R are its own ids and their CHM13 offsets are read off its CHM13 path; in each other arm the
nodes starting at those offsets are its L and R. (Each arm's graph is its sampled.gbz, the graph its reads were
mapped to; `vg convert -f` gives its node lengths and CHM13's walk.)

A read is in a region when its alignment (giraffe's GAF) visits one of the region's nodes in any arm. Over those
reads, per arm: mapped, mean MAPQ, MAPQ 60, mapped at MAPQ 0, mean alignment score and identity; then the arms
pairwise by score, and per region (with at least MIN_READS reads) the arms' MAPQ-0 counts and mean scores, the
regions where they differ most listed. Writes <work dir>/region_reads.tsv (per region and arm).
"""
import bisect
import collections
import gzip
import json
import os
import subprocess
import sys
from array import array

MIN_READS = 20
SHOWN = 15


def node_lengths_and_chm13(vg, gbz):
    """The graph's node lengths ({id: length} as an array from its smallest id) and CHM13's walk [(id, offset)]."""
    proc = subprocess.Popen([vg, "convert", "-f", gbz], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    lengths = {}
    walk = None
    for line in proc.stdout:
        if line.startswith("S\t"):
            _, nid, seq = line.split("\t", 3)[:3]
            lengths[int(nid)] = len(seq.rstrip("\n"))
        elif line.startswith("W\t"):
            w = line.rstrip("\n").split("\t")
            if w[1] == "CHM13" and w[3].startswith("chr20"):
                steps = w[6]
                walk = walk or []
                off = int(w[4])
                k = 0
                while k < len(steps):
                    j = k + 1
                    while j < len(steps) and steps[j] not in "<>":
                        j += 1
                    nid = int(steps[k + 1 : j])
                    walk.append((nid, off))
                    off += lengths[nid]
                    k = j
        elif line.startswith("P\t") and walk is None:
            w = line.rstrip("\n").split("\t")
            if w[1].startswith("CHM13") and "chr20" in w[1]:
                walk = []
                off = 0
                for step in w[2].split(","):
                    nid = int(step[:-1])
                    walk.append((nid, off))
                    off += lengths[nid]
    proc.wait()
    return walk


def region_intervals(packages_path, walk_first, walks):
    """{arm index: sorted [(first id, last id, region)]}: each region's node ids in each arm."""
    offset_of = {nid: off for nid, off in walk_first}
    by_offset = [{off: nid for nid, off in walk} for walk in walks]
    out = [[] for _ in walks]
    with gzip.open(packages_path, "rt") as f:
        for line in f:
            p = json.loads(line)
            L, R = p.get("L"), p.get("R")
            if L is None or R is None or L not in offset_of or R not in offset_of:
                continue
            for a, table in enumerate(by_offset):
                l2, r2 = table.get(offset_of[L]), table.get(offset_of[R])
                if l2 is None or r2 is None:
                    continue
                lo, hi = min(l2, r2), max(l2, r2)
                out[a].append((lo, hi, p["id"]))
    for x in out:
        x.sort()
    return out


def read_regions(gaf, intervals):
    """{read: region} for the reads whose alignment visits a region's nodes (the read keyed as compare_reads.sh
    keys it: name/mate)."""
    starts = [x[0] for x in intervals]
    out = {}
    with gzip.open(gaf, "rt") as f:
        for line in f:
            w = line.split("\t", 13)
            if w[5] == "*":
                continue
            tag = "1" if "\tfn:Z:" in line else "2" if "\tfp:Z:" in line else "s"
            path = w[5]
            k = 0
            hit = None
            while k < len(path) and hit is None:
                j = k + 1
                while j < len(path) and path[j] not in "<>":
                    j += 1
                nid = int(path[k + 1 : j])
                i = bisect.bisect_right(starts, nid) - 1
                if i >= 0 and nid <= intervals[i][1]:
                    hit = intervals[i][2]
                k = j
            if hit is not None:
                out[w[0] + "/" + tag] = hit
    return out


def main():
    work, packages_path, vg = sys.argv[1:4]
    arms = sys.argv[4:]
    walks = [node_lengths_and_chm13(vg, os.path.join(work, "arm_%s" % arm, "sampled.gbz")) for arm in arms]
    intervals = region_intervals(packages_path, walks[0], walks)
    print("regions placed per arm:", ", ".join("%s %d" % (arm, len(x)) for arm, x in zip(arms, intervals)))
    region_of = {}
    for a, arm in enumerate(arms):
        for read, rid in read_regions(os.path.join(work, "arm_%s" % arm, "short.gaf.gz"), intervals[a]).items():
            region_of.setdefault(read, rid)
    print("reads in a region in some arm: %d" % len(region_of))
    k = len(arms)
    stats = [collections.Counter() for _ in arms]
    per_region = collections.defaultdict(lambda: [collections.Counter() for _ in arms])
    versus = collections.Counter()
    with open(os.path.join(work, "reads.joined.tsv")) as f:
        for line in f:
            name = line[: line.index("\t")]
            rid = region_of.get(name)
            if rid is None:
                continue
            w = line.rstrip("\n").split("\t")
            vals = [tuple(int(x) for x in w[1 + 4 * a : 5 + 4 * a]) for a in range(k)]
            for a, (mapq, matches, length, score) in enumerate(vals):
                for s in (stats[a], per_region[rid][a]):
                    s["reads"] += 1
                    s["mapped"] += length > 0
                    s["mapq"] += mapq
                    s["q60"] += mapq >= 60
                    s["q0"] += length > 0 and mapq == 0
                    s["score"] += score
                    s["matches"] += matches
                    s["length"] += length
            for a in range(k):
                for b in range(a + 1, k):
                    sa, sb = vals[a][3], vals[b][3]
                    versus[(a, b, "higher" if sb > sa else "lower" if sb < sa else "same")] += 1
    print("over those reads:")
    for arm, s in zip(arms, stats):
        n = max(s["reads"], 1)
        print("  %-10s reads %d mapped %d, mean MAPQ %.2f, MAPQ 60 %d, mapped at MAPQ 0 %d, mean score %.2f, identity "
              "%.4f" % (arm, s["reads"], s["mapped"], s["mapq"] / n, s["q60"], s["q0"], s["score"] / n,
                        s["matches"] / max(s["length"], 1)))
    for a in range(k):
        for b in range(a + 1, k):
            print("  %s against %s: higher score on %d reads, lower on %d, same on %d" % (
                arms[b], arms[a], versus[(a, b, "higher")], versus[(a, b, "lower")], versus[(a, b, "same")]))
    rows = []
    with open(os.path.join(work, "region_reads.tsv"), "w") as f:
        f.write("region\tarm\treads\tmapped\tmean_mapq\tq60\tq0\tmean_score\tidentity\n")
        for rid, ss in per_region.items():
            for arm, s in zip(arms, ss):
                n = max(s["reads"], 1)
                f.write("%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%.3f\t%.5f\n" % (
                    rid, arm, s["reads"], s["mapped"], s["mapq"] / n, s["q60"], s["q0"], s["score"] / n,
                    s["matches"] / max(s["length"], 1)))
            if ss[0]["reads"] >= MIN_READS:
                rows.append((rid, ss))
    for b in range(1, k):
        print("regions (%d reads or more) where %s has the most more reads at MAPQ 0 than %s, and the most fewer:" % (
            MIN_READS, arms[b], arms[0]))
        rows.sort(key=lambda x: x[1][b]["q0"] - x[1][0]["q0"])
        for rid, ss in rows[:5] + rows[-SHOWN:]:
            print("  %-22s reads %5d | MAPQ 0: %s | mean score: %s" % (
                rid, ss[0]["reads"], " ".join("%s %d" % (arm, s["q0"]) for arm, s in zip(arms, ss)),
                " ".join("%s %.1f" % (arm, s["score"] / max(s["reads"], 1)) for arm, s in zip(arms, ss))))


if __name__ == "__main__":
    main()
