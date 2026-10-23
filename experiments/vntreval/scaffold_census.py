#!/usr/bin/env python3
"""How many columns the haplotypes' optimal pairwise alignments agree on inside the repeat, letting some haplotypes
lack a letter there.

usage: python3 -I experiments/vntreval/scaffold_census.py <vntr-eval> <karan> <regions>

For each region: every pair of distinct haplotypes aligned optimally (karan, unit costs, indels placed left), the
letters they match (=) joined (union-find). A group is a scaffold column when it holds at most one letter of each
haplotype, every two of its members' alignment matches their letters, and every absent haplotype's alignment to
every member leaves the member's letter unmatched. Prints, by the share of haplotypes a column holds, how many such
columns lie in the CHM13 core and outside it, and how many groups fail (two letters of one haplotype, a missing
match, or an absent haplotype matched).
"""
import collections
import json
import os
import re
import subprocess
import sys
import tempfile


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


def main():
    vntr_eval, karan, regions = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
    for rid in regions:
        meta = json.load(open(os.path.join(vntr_eval, "regions", rid, "region.json")))
        names, seqs = read_fasta(os.path.join(vntr_eval, "regions", rid, "hap32.fa"))
        chm13 = next(k for k, nm in enumerate(names) if nm.upper().startswith("CHM13"))
        core = (meta["core_start"] - meta["span_start"], meta["core_end"] - meta["span_start"])
        order = list(dict.fromkeys(seqs))
        ref = order.index(seqs[chm13])
        n = len(order)
        pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
        with tempfile.TemporaryDirectory() as tmp:
            tsv = os.path.join(tmp, "p.tsv")
            with open(tsv, "w") as f:
                for k, (a, b) in enumerate(pairs):
                    f.write("%d\t%s\t%s\n" % (k, order[b], order[a]))
            out = subprocess.run([karan, "--pairs", tsv], capture_output=True, text=True).stdout
        maps = {}
        for line in out.strip().splitlines()[1:]:
            f = line.split("\t")
            a, b = pairs[int(f[0])]
            m = [-1] * len(order[a])
            i = j = 0
            for k, op in re.findall(r"(\d+)([=XID])", f[9]):
                for _ in range(int(k)):
                    if op == "D":
                        i += 1
                    elif op == "I":
                        j += 1
                    else:
                        if op == "=":
                            m[j] = i
                        i += 1
                        j += 1
            maps[(a, b)] = m
        offset = [0]
        for s in order:
            offset.append(offset[-1] + len(s))
        parent = list(range(offset[-1]))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for (a, b), m in maps.items():
            if a < b:
                for i, j in enumerate(m):
                    if j >= 0:
                        x, y = find(offset[a] + i), find(offset[b] + j)
                        if x != y:
                            parent[x] = y
        groups = collections.defaultdict(list)
        for s in range(n):
            for i in range(len(order[s])):
                groups[find(offset[s] + i)].append((s, i))
        census = collections.Counter()
        for members in groups.values():
            held = {s for s, _ in members}
            if len(members) < 2:
                continue
            if len(held) != len(members):
                census["fail: two letters of one haplotype"] += 1
                continue
            letter = dict(members)
            ok = all(maps[(a, b)][letter[a]] == letter[b] for a in held for b in held if a != b)
            if not ok:
                census["fail: a member pair does not match"] += 1
                continue
            absent_ok = all(maps[(a, x)][letter[a]] == -1 for a in held for x in range(n) if x not in held)
            if not absent_ok:
                census["fail: an absent haplotype is matched"] += 1
                continue
            share = len(held) / n
            band = "all" if share == 1 else ">=0.9" if share >= 0.9 else ">=0.5" if share >= 0.5 else "<0.5"
            where = "flank"
            if ref in letter and core[0] <= letter[ref] < core[1]:
                where = "core"
            elif ref not in letter:
                where = "not in CHM13"
            census["%s, %s" % (band, where)] += 1
        print("## %s (%d distinct haplotypes, core %d bp)" % (rid, n, core[1] - core[0]))
        for k in sorted(census):
            print("  %-45s %d" % (k, census[k]))
        # Star columns: from the medoid's letters, members pruned until every pair agrees.
        cost = [sum(sum(1 for x in maps[(a, b)] if x < 0) for b in range(n) if b != a) for a in range(n)]
        centre = min(range(n), key=lambda a: cost[a])
        star = collections.Counter()
        for j in range(len(order[centre])):
            letter = {centre: j}
            for x in range(n):
                if x != centre and maps[(centre, x)][j] >= 0:
                    letter[x] = maps[(centre, x)][j]
            while len(letter) > 1:
                agree = {a: sum(maps[(a, b)][letter[a]] == letter[b] for b in letter if b != a) for a in letter}
                worst = min(agree, key=lambda a: (agree[a], a != centre))
                if agree[worst] == len(letter) - 1:
                    break
                del letter[worst]
            share = len(letter) / n
            band = "all" if share == 1 else ">=0.9" if share >= 0.9 else ">=0.5" if share >= 0.5 else "<0.5"
            where = "flank"
            if ref in letter and core[0] <= letter[ref] < core[1]:
                where = "core"
            elif ref not in letter:
                where = "not in CHM13"
            star["%s, %s" % (band, where)] += 1
        print("  star columns from the medoid (haplotype %d):" % centre)
        for k in sorted(star):
            print("    %-43s %d" % (k, star[k]))


if __name__ == "__main__":
    main()
