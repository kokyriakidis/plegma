#!/usr/bin/env python3
"""The reads whose alignments differ between chr20 arms (compare_reads.sh's joined table: per read, per arm,
mapping quality, matches, aligned length and alignment score), and which arm aligns them better.

usage: python3 -I experiments/chr20/compare_reads.py <reads.joined.tsv> <arm> <arm> [<arm> ...]

Prints how many reads differ in any arm; per arm, over those reads, how many it maps, its mean alignment score,
identity (matches over aligned length) and mapping quality, how many it maps at MAPQ 60 and at MAPQ 0; and for each
pair of arms, on how many differing reads the second scores higher, lower or the same.
"""
import itertools
import sys


def main():
    path = sys.argv[1]
    arms = sys.argv[2:]
    k = len(arms)
    total = differ = 0
    stats = [dict(mapped=0, score=0, matches=0, length=0, mapq=0, q60=0, q0=0) for _ in arms]
    versus = {pair: [0, 0, 0] for pair in itertools.combinations(range(k), 2)}
    with open(path) as f:
        for line in f:
            w = line.rstrip("\n").split("\t")
            total += 1
            vals = [tuple(int(x) for x in w[1 + 4 * a : 5 + 4 * a]) for a in range(k)]
            if all(v == vals[0] for v in vals[1:]):
                continue
            differ += 1
            for a, (mapq, matches, length, score) in enumerate(vals):
                s = stats[a]
                if length > 0:
                    s["mapped"] += 1
                s["score"] += score
                s["matches"] += matches
                s["length"] += length
                s["mapq"] += mapq
                s["q60"] += mapq >= 60
                s["q0"] += length > 0 and mapq == 0
            for (a, b), t in versus.items():
                sa, sb = vals[a][3], vals[b][3]
                t[0 if sb > sa else 1 if sb < sa else 2] += 1
    print("reads %d, differing in some arm %d (%.3f%%)" % (total, differ, 100.0 * differ / max(total, 1)))
    print("over the differing reads:")
    for arm, s in zip(arms, stats):
        print("  %-10s mapped %d, mean score %.2f, identity %.4f, mean MAPQ %.2f, MAPQ 60 %d, mapped at MAPQ 0 %d" % (
            arm, s["mapped"], s["score"] / max(differ, 1), s["matches"] / max(s["length"], 1),
            s["mapq"] / max(differ, 1), s["q60"], s["q0"]))
    for (a, b), (better, worse, same) in versus.items():
        print("  %s against %s: higher score on %d reads, lower on %d, same on %d" % (arms[b], arms[a], better, worse,
                                                                                    same))


if __name__ == "__main__":
    main()
