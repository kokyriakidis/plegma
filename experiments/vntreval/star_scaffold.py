#!/usr/bin/env python3
"""A prototype starting MSA from the medoid's star columns, pruned to agreement, for plegma's refinement (-i, -I).

usage: python3 -I experiments/vntreval/star_scaffold.py <karan> <in.fa> <out.msa.fa> [min share]

The medoid (the distinct sequence nearest all others in summed edit distance) is aligned optimally to every other
(karan, unit costs, indels placed left). For each of its letters, the letters the alignments match it to form a
candidate column, pruned (the member agreeing with fewest others first) until every two members' own optimal
alignment matches their letters; a column keeping at least `min share` (0.5) of the sequences is a scaffold column.
A sequence's letters outside scaffold columns go between them where its alignment to the medoid puts them, each
sequence's own run of them in columns of its own: the windows between scaffold columns are left for the
refinement to align. Identical sequences share their row.
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


def main():
    karan, src, dst = sys.argv[1:4]
    share = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
    names, seqs = read_fasta(src)
    order = list(dict.fromkeys(seqs))
    n = len(order)
    pairs = [(a, b) for a in range(n) for b in range(n) if a != b]
    with tempfile.TemporaryDirectory() as tmp:
        tsv = os.path.join(tmp, "p.tsv")
        with open(tsv, "w") as f:
            for k, (a, b) in enumerate(pairs):
                f.write("%d\t%s\t%s\n" % (k, order[b], order[a]))
        out = subprocess.run([karan, "--pairs", tsv], capture_output=True, text=True).stdout
    maps, costs = {}, {}
    for line in out.strip().splitlines()[1:]:
        f = line.split("\t")
        a, b = pairs[int(f[0])]
        costs[(a, b)] = int(f[3])
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
    centre = min(range(n), key=lambda a: sum(costs[(a, b)] for b in range(n) if b != a))
    # Scaffold columns: centre letter j -> {sequence: letter}.
    scaffold = []
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
        if len(letter) >= share * n:
            scaffold.append(letter)
    # Each sequence's letters between consecutive scaffold columns it is in, placed by its alignment to the centre:
    # a letter goes in the window after the last scaffold column whose centre letter precedes it.
    # windows[w][s]: sequence s's letters in window w (before scaffold column w; w = len(scaffold) after the last).
    windows = [[[] for _ in range(n)] for _ in range(len(scaffold) + 1)]
    import bisect
    centre_of = [letter[centre] for letter in scaffold]
    for s in range(n):
        mine = [(letter[s], w) for w, letter in enumerate(scaffold) if s in letter]
        mine_letters = [i for i, _ in mine]
        m = maps[(s, centre)] if s != centre else list(range(len(order[centre])))
        last_centre = -1
        for i in range(len(order[s])):
            if m[i] >= 0:
                last_centre = m[i]
            k = bisect.bisect_left(mine_letters, i)
            if k < len(mine) and mine_letters[k] == i:
                continue
            # Allowed windows: after s's scaffold column before i, up to its scaffold column after i.
            lo = mine[k - 1][1] + 1 if k > 0 else 0
            hi = mine[k][1] if k < len(mine) else len(scaffold)
            # Where the centre puts it: after every scaffold column whose centre letter is at or before the one i faces.
            by_centre = bisect.bisect_right(centre_of, last_centre)
            windows[max(lo, min(hi, by_centre))][s].append(order[s][i])
    rows = [[] for _ in range(n)]
    for w in range(len(scaffold) + 1):
        width = max(len(p) for p in windows[w][:]) if n else 0
        for s in range(n):
            piece = "".join(windows[w][s])
            rows[s].append(piece + "-" * (width - len(piece)))
        if w < len(scaffold):
            for s in range(n):
                rows[s].append(order[s][scaffold[w][s]] if s in scaffold[w] else "-")
    row_of = {order[s]: "".join(rows[s]) for s in range(n)}
    for k, s in enumerate(seqs):
        assert row_of[s].replace("-", "") == s, "row %d does not spell its sequence" % k
    with open(dst, "w") as f:
        for name, s in zip(names, seqs):
            f.write(">%s\n%s\n" % (name, row_of[s]))
    print("centre %d, scaffold columns %d of its %d letters" % (centre, len(scaffold), len(order[centre])),
          file=sys.stderr)


if __name__ == "__main__":
    main()
