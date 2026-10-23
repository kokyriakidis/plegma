#!/usr/bin/env python3
"""Whether any MSA could pair a pair's copies as its optimal alignment does: the pair aligned through each third
sequence x (both optimally to x, merged through x's letters, the three-row MSA an MSA aligning both well to x
implies), and that MSA's copy pairing of the pair audited by build/structure.

usage: python3 -I experiments/anchors/star_copies.py <structure> <karan> <sequences.fa> <a> <b> <period> <out dir>
           [--affine]

Prints, for each x, the share of a's copies the three-row MSA pairs right (structure's `right`), against the region
MSA's share.
"""
import os
import re
import subprocess
import sys


def read(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return list(dict.fromkeys("".join(s) for s in seqs))


def pairing(karan, ref, query, work):
    """Each letter of ref: its partner in query under karan's optimal unit-cost alignment (-1: deleted)."""
    tsv = os.path.join(work, "pair.tsv")
    with open(tsv, "w") as f:
        f.write("q\t%s\t%s\n" % (ref, query))
    out = subprocess.run([karan, "--pairs", tsv, "--mode", "global"], capture_output=True, text=True).stdout
    line = out.strip().splitlines()[-1].split("\t")
    i, j = int(line[5]), int(line[7])
    partner = [-1] * len(ref)
    for n, op in re.findall(r"(\d+)([=XIDM])", line[9]):
        for _ in range(int(n)):
            if op == "D":
                i += 1
            elif op == "I":
                j += 1
            else:
                partner[i] = j
                i += 1
                j += 1
    return partner


def star(x, a, b, to_a, to_b):
    """Three rows: x, a, b, merged through x's letters (a's and b's letters x does not pair in columns of their own,
    in order)."""
    rows = {"x": [], "a": [], "b": []}
    pa = pb = 0
    for i in range(len(x) + 1):
        ja = to_a[i] if i < len(x) else len(a)
        jb = to_b[i] if i < len(x) else len(b)
        # Letters of a and b before this column, unpaired: own columns.
        if ja >= 0:
            while pa < ja:
                rows["x"].append("-"); rows["a"].append(a[pa]); rows["b"].append("-"); pa += 1
        if jb >= 0:
            while pb < jb:
                rows["x"].append("-"); rows["a"].append("-"); rows["b"].append(b[pb]); pb += 1
        if i == len(x):
            break
        rows["x"].append(x[i])
        if ja >= 0:
            rows["a"].append(a[ja]); pa = ja + 1
        else:
            rows["a"].append("-")
        if jb >= 0:
            rows["b"].append(b[jb]); pb = jb + 1
        else:
            rows["b"].append("-")
    return {k: "".join(v) for k, v in rows.items()}


def right_share(structure, fa, msa, period, affine):
    args = [structure, fa, msa, str(period), "--pair", "1", "2"] + (["--affine"] if affine else [])
    out = subprocess.run(args, capture_output=True, text=True).stdout
    m = re.search(r"right: (\d+) \(([\d.]+)%\)", out)
    return float(m.group(2)) if m else None


def main():
    args = sys.argv[1:]
    affine = "--affine" in args
    args = [x for x in args if x != "--affine"]
    structure, karan, fa, ia, ib, period, out = args
    ia, ib, period = int(ia), int(ib), int(period)
    os.makedirs(out, exist_ok=True)
    seqs = read(fa)
    a, b = seqs[ia], seqs[ib]
    shares = []
    for x in range(len(seqs)):
        if x in (ia, ib):
            continue
        c = seqs[x]
        rows = star(c, a, b, pairing(karan, c, a, out), pairing(karan, c, b, out))
        three = os.path.join(out, "three.%d.fa" % x)
        with open(three, "w") as f:
            f.write(">x\n%s\n>a\n%s\n>b\n%s\n" % (c, a, b))
        msa = os.path.join(out, "three.%d.msa.fa" % x)
        with open(msa, "w") as f:
            for k in ("x", "a", "b"):
                f.write(">%s\n%s\n" % (k, rows[k]))
        # structure's pair: a and b are sequences 1 and 2 of the three (x is 0), the longer of them taken as a.
        share = right_share(structure, three, msa, period, affine)
        shares.append((share, x))
        print("through %2d: copies right %s%%" % (x, share), flush=True)
    shares = [s for s in shares if s[0] is not None]
    shares.sort(reverse=True)
    print("best %s%% (through %d), median %s%%" % (shares[0][0], shares[0][1], shares[len(shares) // 2][0]))


if __name__ == "__main__":
    main()
