#!/usr/bin/env python3
"""Where two MSAs of one region differ, by how much: whole copies of the repeat unit, or less.

usage: python3 -I experiments/vntreval/copy_shift.py <msa X> <msa Y> <period>

For every pair of distinct sequences (a, b) and every letter of a, the letter of b it faces in X and in Y (or a
gap). Letters facing the same partner in both are `same`; letters facing a letter in both but different ones are
classed by the offset d between the two partners: `copies k` when d is within 2 letters of k periods (k != 0),
`near` when |d| <= 2, `other` otherwise; `gap/letter` when one MSA gaps the partner and the other does not. Prints the
shares of each class over all pairs, the commonest copy offsets, and the pairs differing most.
"""
import collections
import sys


def read_rows(path):
    names, rows = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            rows.append([])
        elif line:
            rows[-1].append(line.upper().replace(".", "-"))
    return names, ["".join(r) for r in rows]


def by_sequence(path):
    out = {}
    for n, r in zip(*read_rows(path)):
        out.setdefault(r.replace("-", ""), (n, r))
    return out


def facing(ra, rb):
    """For each letter of a, the index of b's letter in its column, or -1."""
    out = []
    j = 0
    for x, y in zip(ra, rb):
        if x != "-":
            out.append(j if y != "-" else -1)
        if y != "-":
            j += 1
    return out


def main():
    x, y, p = by_sequence(sys.argv[1]), by_sequence(sys.argv[2]), int(sys.argv[3])
    seqs = [s for s in x if s in y]
    total = collections.Counter()
    shifts = collections.Counter()
    worst = []
    for i in range(len(seqs)):
        for k in range(len(seqs)):
            if i == k:
                continue
            fx = facing(x[seqs[i]][1], x[seqs[k]][1])
            fy = facing(y[seqs[i]][1], y[seqs[k]][1])
            pair = collections.Counter()
            for a, b in zip(fx, fy):
                if a == b:
                    pair["same"] += 1
                elif a < 0 or b < 0:
                    pair["gap/letter"] += 1
                else:
                    d = a - b
                    c = round(d / p)
                    if c != 0 and abs(d - c * p) <= 2:
                        pair["copies"] += 1
                        shifts[c] += 1
                    elif abs(d) <= 2:
                        pair["near"] += 1
                    else:
                        pair["other"] += 1
            total.update(pair)
            worst.append((sum(v for key, v in pair.items() if key != "same"), x[seqs[i]][0], x[seqs[k]][0], pair))
    n = sum(total.values())
    print("%d sequences, %d letter-pairs; period %d" % (len(seqs), n, p))
    for key in ("same", "copies", "near", "other", "gap/letter"):
        print("  %-11s %6.2f%%" % (key, 100.0 * total[key] / max(1, n)))
    print("copy offsets (X minus Y, in copies):", dict(shifts.most_common(8)))
    worst.sort(reverse=True)
    for w, a, b, pair in worst[:6]:
        print("  %-28s %-28s %s" % (a[:28], b[:28], dict(pair)))


if __name__ == "__main__":
    main()
