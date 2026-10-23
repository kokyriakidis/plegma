#!/usr/bin/env python3
"""One region, one pair: the flanks pinned first, then the seeds inside the array, against the optimum and an MSA.

usage: python3 -I experiments/vntreval/one_region.py <karan> <hap32.fa> <region.json> <name a> <name b>
           [msa.fa] [k] [step]

1. Flanks: k-mers (21 by default) unique in each whole sequence and outside tandem repeats (a 12-mer recurring
   within 500 letters masks it), shared by the two, chained along rising positions: the last before the array
   and the first after it pin the array window.
2. Inside the window: k-mers unique within each sequence's array piece and shared, chained (longest chain rising in
   both): the copies a seed sits in carry something of their own, so the pair of copies is homologous.
3. Along a, every `step` letters (250): the diagonal (b's offset minus a's) of the window seeds' chain, of the
   optimal unit-cost alignment (karan), and of the MSA's implied alignment when an MSA is given.
"""
import bisect
import json
import re
import subprocess
import sys
import tempfile


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line.upper())
    return dict(zip(names, ["".join(s) for s in seqs]))


def kmer_positions(s, k):
    pos = {}
    for i in range(len(s) - k + 1):
        pos.setdefault(s[i:i + k], []).append(i)
    return pos


def masked(s, k=12, distance=500):
    """Positions whose 12-mer recurs within `distance` letters."""
    last, out = {}, [False] * len(s)
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if w in last and i - last[w] <= distance:
            out[i] = out[last[w]] = True
        last[w] = i
    return out


def chain(seeds):
    """The longest chain of (i, j) seeds rising in both."""
    seeds.sort()
    tails, back, idx = [], [-1] * len(seeds), []
    for t, (_, j) in enumerate(seeds):
        k = bisect.bisect_left([seeds[x][1] for x in idx], j)
        if k > 0:
            back[t] = idx[k - 1]
        if k == len(idx):
            idx.append(t)
        else:
            idx[k] = t
    out, t = [], idx[-1] if idx else -1
    while t >= 0:
        out.append(seeds[t])
        t = back[t]
    return out[::-1]


def unique_shared(a, b, k, a_ok=None, b_ok=None):
    pa, pb = kmer_positions(a, k), kmer_positions(b, k)
    seeds = []
    for w, ia in pa.items():
        if len(ia) == 1 and len(pb.get(w, ())) == 1:
            i, j = ia[0], pb[w][0]
            if (a_ok is None or a_ok(i)) and (b_ok is None or b_ok(j)):
                seeds.append((i, j))
    return seeds


def optimal_pairing(karan, a, b):
    with tempfile.TemporaryDirectory() as tmp:
        open(tmp + "/a.fa", "w").write(">a\n%s\n" % a)
        open(tmp + "/b.fa", "w").write(">b\n%s\n" % b)
        cigar = subprocess.run([karan, tmp + "/a.fa", tmp + "/b.fa"], capture_output=True,
                               text=True).stdout.strip().splitlines()[-1].split("\t")[9]
    pair, i, j = [-1] * len(a), 0, 0
    for n, op in re.findall(r"(\d+)([=XID])", cigar):
        for _ in range(int(n)):
            if op == "D":
                i += 1
            elif op == "I":
                j += 1
            else:
                pair[i] = j
                i += 1
                j += 1
    return pair


def msa_pairing(path, na, nb, n):
    rows = read_fasta(path)
    ra, rb = rows[na], rows[nb]
    pair, i, j = [-1] * n, 0, 0
    for x, y in zip(ra, rb):
        if x != "-" and y != "-":
            pair[i] = j
        i += x != "-"
        j += y != "-"
    return pair


def main():
    karan, fasta, meta_path, na, nb = sys.argv[1:6]
    msa = sys.argv[6] if len(sys.argv) > 6 and sys.argv[6] != "-" else None
    k = int(sys.argv[7]) if len(sys.argv) > 7 else 21
    step = int(sys.argv[8]) if len(sys.argv) > 8 else 250
    seqs = read_fasta(fasta)
    a, b = seqs[na], seqs[nb]
    meta = json.load(open(meta_path))
    period = meta.get("period") or 0
    # 1. Flanks.
    ma, mb = masked(a), masked(b)
    flank = chain(unique_shared(a, b, k, lambda i: not any(ma[i:i + k]), lambda j: not any(mb[j:j + k])))
    print("1. flank seeds (unique in the whole sequences, outside repeats): %d" % len(flank))
    if not flank:
        print("   none: no window")
        return
    # The array window: the widest gap between consecutive flank seeds.
    gaps = [(flank[t + 1][0] - flank[t][0], t) for t in range(len(flank) - 1)]
    _, t = max(gaps)
    (la, lb), (ra_, rb_) = flank[t], flank[t + 1]
    wa0, wa1, wb0, wb1 = la + k, ra_, lb + k, rb_
    print("   array window: a %d..%d (%d letters), b %d..%d (%d letters); %d flank seeds left of it, %d right" % (
        wa0, wa1, wa1 - wa0, wb0, wb1, wb1 - wb0, t + 1, len(flank) - t - 1))
    # 2. Inside the window.
    pa_, pb_ = a[wa0:wa1], b[wb0:wb1]
    inside = chain([(i + wa0, j + wb0) for i, j in unique_shared(pa_, pb_, k)])
    print("2. window seeds (unique within the two array pieces), chained: %d" % len(inside))
    # 3. Diagonals along a.
    opt = optimal_pairing(karan, a, b)
    imp = msa_pairing(msa, na, nb, len(a)) if msa else None
    print("3. %8s %10s %10s %10s %s" % ("a pos", "seeds", "optimum", "MSA" if msa else "", "(diagonal: b - a; units of %s)" % period))
    seed_i = [i for i, _ in inside]

    def diag_seed(p):
        x = bisect.bisect_right(seed_i, p) - 1
        return inside[x][1] - inside[x][0] if x >= 0 else None

    def diag_pair(pair, p):
        while p >= 0 and pair[p] < 0:
            p -= 1
        return pair[p] - p if p >= 0 else None

    def show(d):
        if d is None:
            return "-"
        return "%d (%+.1f)" % (d, d / period) if period else str(d)

    for p in range(wa0 - step, wa1 + step, step):
        if p < 0 or p >= len(a):
            continue
        tag = "" if wa0 <= p < wa1 else " flank"
        print("   %8d %14s %14s %14s%s" % (p, show(diag_seed(p)), show(diag_pair(opt, p)),
                                       show(diag_pair(imp, p)) if imp else "", tag))


if __name__ == "__main__":
    main()
