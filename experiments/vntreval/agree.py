#!/usr/bin/env python3
"""A prototype MSA built from where the haplotypes' optimal pairwise alignments agree, recursing into the rest.

usage: python3 -I experiments/vntreval/agree.py <karan> <plegma> <in.fa> <out.msa.fa> [plegma flags for leaves...]

1. Every pair of distinct sequences is aligned optimally (karan, unit costs, indels placed left).
2. Letters their alignments pair (=) are joined (union-find). A group holding exactly one letter of every sequence
   is a settled column; the settled columns kept are those whose letters rise in every sequence together (greedy,
   in the first sequence's order).
3. The settled columns split the sequences into windows: each is solved again from step 1, on its pieces alone
   (their own optimal alignments, with both ends fixed, can agree where the whole sequences' did not).
4. A window where no column settles is a leaf: its pieces are aligned by `plegma [flags] -r 1` (`-L -I 20,64` by
   default).
Identical sequences are aligned once and share their row. Prints the windows' depths and sizes to stderr.
"""
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


class Solver:
    def __init__(self, karan, plegma, flags, tmp):
        self.karan, self.plegma, self.flags, self.tmp = karan, plegma, flags, tmp
        self.leaves = []
        self.settled = 0

    def pairwise(self, seqs):
        """maps[(a, b)][i]: the letter of b that letter i of a is matched (=) to, or -1."""
        n = len(seqs)
        pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
        tsv = os.path.join(self.tmp, "pairs.tsv")
        with open(tsv, "w") as f:
            for k, (a, b) in enumerate(pairs):
                f.write("%d\t%s\t%s\n" % (k, seqs[b] or "N", seqs[a] or "N"))
        out = subprocess.run([self.karan, "--pairs", tsv], capture_output=True, text=True).stdout
        maps = {}
        for line in out.strip().splitlines()[1:]:
            f = line.split("\t")
            a, b = pairs[int(f[0])]
            m = [-1] * len(seqs[a])
            i = j = 0
            for k, op in re.findall(r"(\d+)([=XID])", f[9]):
                for _ in range(int(k)):
                    if op == "D":
                        i += 1
                    elif op == "I":
                        j += 1
                    else:
                        if op == "=" and seqs[a] and seqs[b]:
                            m[j] = i
                        i += 1
                        j += 1
            maps[(a, b)] = m
        return maps

    def settled_columns(self, seqs):
        n = len(seqs)
        offset = [0]
        for s in seqs:
            offset.append(offset[-1] + len(s))
        parent = list(range(offset[-1]))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        maps = self.pairwise(seqs)
        for (a, b), m in maps.items():
            for i, j in enumerate(m):
                if j >= 0:
                    x, y = find(offset[a] + i), find(offset[b] + j)
                    if x != y:
                        parent[x] = y
        groups = {}
        for s in range(n):
            for i in range(len(seqs[s])):
                groups.setdefault(find(offset[s] + i), []).append((s, i))
        columns = []
        for members in groups.values():
            if len(members) == n and len({s for s, _ in members}) == n:
                pos = [0] * n
                for s, i in members:
                    pos[s] = i
                columns.append(pos)
        columns.sort()
        kept = []
        last = [-1] * n
        for pos in columns:
            if all(p > q for p, q in zip(pos, last)):
                kept.append(pos)
                last = pos
        return kept

    def leaf(self, seqs):
        n = len(seqs)
        longest = max(len(s) for s in seqs)
        if longest == 0:
            return [""] * n
        if sum(1 for s in seqs if s) == 1 or len(set(seqs)) == 1:
            return [s.ljust(longest, "-") if s else "-" * longest for s in seqs]
        path = os.path.join(self.tmp, "leaf.fa")
        with open(path, "w") as f:
            for k, s in enumerate(seqs):
                f.write(">s%d\n%s\n" % (k, s or ""))
        out = subprocess.run([self.plegma] + self.flags + ["-r", "1", path], capture_output=True, text=True).stdout
        rows = {}
        name = None
        for line in out.splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                rows[name] = []
            elif name:
                rows[name].append(line.strip())
        width = max(len("".join(r)) for r in rows.values()) if rows else 0
        result = []
        for k in range(n):
            r = "".join(rows.get("s%d" % k, []))
            result.append(r if r else "-" * width)
        self.leaves.append((n, longest))
        return result

    def solve(self, seqs, depth=0):
        """Rows of the MSA of `seqs` (all distinct or not)."""
        if all(len(s) == 0 for s in seqs):
            return [""] * len(seqs)
        if any(len(s) == 0 for s in seqs) or depth > 12:
            return self.leaf(seqs)
        columns = self.settled_columns(seqs)
        if not columns:
            return self.leaf(seqs)
        self.settled += len(columns)
        rows = [[] for _ in seqs]
        start = [0] * len(seqs)
        for pos in columns + [None]:
            end = [len(s) for s in seqs] if pos is None else pos
            pieces = [s[a:b] for s, a, b in zip(seqs, start, end)]
            if any(pieces):
                sub = self.solve(pieces, depth + 1)
                for r, x in zip(rows, sub):
                    r.append(x)
            if pos is not None:
                for r, s, p in zip(rows, seqs, pos):
                    r.append(s[p])
                start = [p + 1 for p in pos]
        return ["".join(r) for r in rows]


def main():
    karan, plegma, src, dst = sys.argv[1:5]
    flags = sys.argv[5:] or ["-L", "-I", "20,64"]
    names, seqs = read_fasta(src)
    distinct = list(dict.fromkeys(seqs))
    with tempfile.TemporaryDirectory() as tmp:
        solver = Solver(karan, plegma, flags, tmp)
        rows = solver.solve(distinct)
    row_of = dict(zip(distinct, rows))
    with open(dst, "w") as f:
        for name, s in zip(names, seqs):
            f.write(">%s\n%s\n" % (name, row_of[s]))
    print("settled columns %d, leaves %d, largest leaf %s" % (
        solver.settled, len(solver.leaves), max(solver.leaves, key=lambda x: x[1]) if solver.leaves else None),
        file=sys.stderr)


if __name__ == "__main__":
    main()
