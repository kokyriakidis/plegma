#!/usr/bin/env python3
"""Pairs of VNTR haplotypes with their true alignment.

usage: python3 -I experiments/seeding/simulate_vntr.py <out prefix> <pairs> [seed]

An ancestor is a random flank, `copies` copies of a random unit of `period` letters (each copy carrying its
own variants at `unit_variation` per letter, so copies are similar but not identical) and another flank. Each of
the two haplotypes is derived from it independently: whole copies lost (a run of them) or gained (a run copied
in place, the new letters homologous to nothing in the other haplotype), then substitutions, insertions and
deletions at `mutation` per letter. Every letter keeps its ancestral origin, so the true alignment pairs the
letters that share one.

Writes <prefix>.fa (pair k: records pk_a then pk_b) and <prefix>.truth.tsv: name, period, copies,
unit_variation, mutation, copy-number difference, and the true alignment as a CIGAR of b against a
(=, X, I: a letter of b alone, D: a letter of a alone).
"""
import random
import sys


def mutate(rng, letters, rate, hap):
    out = []
    counter = 0
    for base, key in letters:
        r = rng.random()
        if r < rate / 3:
            out.append((rng.choice([c for c in "ACGT" if c != base]), key))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append((base, key))
            counter += 1
            out.append((rng.choice("ACGT"), (key[0], key[1] + 0.5, hap, counter)))
        else:
            out.append((base, key))
    return out


def derive(rng, flank_left, copies, flank_right, change, rate, hap):
    """One haplotype: `change` copies lost (negative) or gained (positive) in one run, then point mutations."""
    copies = list(copies)
    if change < 0:
        start = rng.randint(0, len(copies) + change)
        del copies[start:start - change]
    elif change > 0:
        start = rng.randint(0, len(copies) - 1)
        run = [copies[(start + k) % len(copies)] for k in range(change)]
        fresh = []
        for k, copy in enumerate(run):
            # Gained letters sort right after the copy they follow and belong to this haplotype alone.
            anchor_key = copies[start + len(run) - 1 if start + len(run) - 1 < len(copies) else -1][-1][1]
            fresh.append([(b, (anchor_key[0], anchor_key[1] + 0.25, hap, k * 10000 + i)) for i, (b, _) in enumerate(copy)])
        copies = copies[:start + len(run)] + fresh + copies[start + len(run):]
    letters = list(flank_left) + [x for copy in copies for x in copy] + list(flank_right)
    return mutate(rng, letters, rate, hap)


def cigar_of(a, b):
    """The alignment pairing letters of equal origin, in origin order."""
    ka = [key for _, key in a]
    kb = [key for _, key in b]
    ops = []
    i = j = 0
    while i < len(a) or j < len(b):
        if i < len(a) and j < len(b) and ka[i] == kb[j]:
            ops.append("=" if a[i][0] == b[j][0] else "X")
            i += 1
            j += 1
        elif j >= len(b) or (i < len(a) and ka[i] < kb[j]):
            ops.append("D")
            i += 1
        else:
            ops.append("I")
            j += 1
    out = []
    run = 1
    for k in range(1, len(ops) + 1):
        if k < len(ops) and ops[k] == ops[k - 1]:
            run += 1
        else:
            out.append("%d%s" % (run, ops[k - 1]))
            run = 1
    return "".join(out)


def main():
    prefix, pairs = sys.argv[1], int(sys.argv[2])
    rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
    with open(prefix + ".fa", "w") as fa, open(prefix + ".truth.tsv", "w") as tsv:
        tsv.write("name\tperiod\tcopies\tunit_variation\tmutation\tcopy_difference\tcigar\n")
        for k in range(pairs):
            period = rng.choice([2, 3, 6, 9, 18, 30, 34, 45, 60])
            copies_n = rng.randint(5, 2000 // period + 5)
            unit_variation = rng.choice([0.0, 0.01, 0.03, 0.08])
            mutation = rng.choice([0.001, 0.005, 0.02])
            pos = 0
            flank_left = []
            for _ in range(rng.randint(150, 400)):
                flank_left.append((rng.choice("ACGT"), (pos, 0.0, 0, 0)))
                pos += 1
            unit = "".join(rng.choice("ACGT") for _ in range(period))
            copies = []
            for _ in range(copies_n):
                copy = []
                for c in unit:
                    base = rng.choice([x for x in "ACGT" if x != c]) if rng.random() < unit_variation else c
                    copy.append((base, (pos, 0.0, 0, 0)))
                    pos += 1
                copies.append(copy)
            flank_right = []
            for _ in range(rng.randint(150, 400)):
                flank_right.append((rng.choice("ACGT"), (pos, 0.0, 0, 0)))
                pos += 1
            spread = max(1, copies_n // 3)
            change_a = rng.randint(-spread, spread) if rng.random() < 0.7 else 0
            change_b = rng.randint(-spread, spread)
            a = derive(rng, flank_left, copies, flank_right, change_a, mutation, 1)
            b = derive(rng, flank_left, copies, flank_right, change_b, mutation, 2)
            name = "p%d" % k
            fa.write(">%s_a\n%s\n>%s_b\n%s\n" % (name, "".join(x for x, _ in a), name, "".join(x for x, _ in b)))
            tsv.write("%s\t%d\t%d\t%.2f\t%.3f\t%d\t%s\n" % (name, period, copies_n, unit_variation, mutation,
                                                          change_b - change_a, cigar_of(a, b)))


if __name__ == "__main__":
    main()
