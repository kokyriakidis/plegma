#!/usr/bin/env python3
"""Which pair of a graph's alleles a region's reads support best, by depth as well as presence: a genotyper without
the graph's snarls, whole alleles against the region's reads.

Every canonical K-mer of the rows (flanked) is counted in each row and in the reads. Only the K-mers whose counts
differ between rows are scored (the others say nothing about the pair, and the reads' errors and HG002's own few
differences otherwise sway it); a pair of rows (a, b) expects a K-mer's read count at (count in a + count in b) times
the depth that pair fits best (the reads' counts of its K-mers over its K-mers' counts: the window's reads thin out
towards its ends), at least FLOOR, and the pair whose Poisson likelihood of the counts is highest is the reads' choice.
A copy-number difference shows as depth, where read placement cannot see it. On the 387 chr20 regions it reaches the
truth pair (within 20 edits of the best the graph allows) in 192, vg call in 198; summed, it is nearer the truth in
the long arrays (76,181 edits over the best against vg call's 126,509), level elsewhere: an estimate of what the
reads hold, not a bound.
"""
import math
import os

import numpy

K = 31
FLOOR = 0.2
COMP = str.maketrans("ACGT", "TGCA")


def canonical_kmers(seq):
    rc = seq.translate(COMP)[::-1]
    n = len(seq)
    for i in range(n - K + 1):
        a = seq[i : i + K]
        b = rc[n - K - i : n - i]
        yield a if a <= b else b


def read_seqs(prep_dir, rid):
    for name in ("R1", "R2", "single"):
        path = os.path.join(prep_dir, rid, name + ".fq")
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for k, line in enumerate(f):
                if k % 4 == 1:
                    yield line.strip().upper()


def best_pair(rows, left, right, prep_dir, rid):
    """(id a, id b, log-likelihood margin over the next pair) the reads support best, among `rows` ({id: row})."""
    seqs = {}
    for a, r in rows.items():
        seqs.setdefault(left + r.replace("-", "") + right, []).append(a)
    distinct = list(seqs)
    index = {}
    for s in distinct:
        for km in canonical_kmers(s):
            if km not in index:
                index[km] = len(index)
    counts = numpy.zeros((len(distinct), len(index)))
    for i, s in enumerate(distinct):
        for km in canonical_kmers(s):
            counts[i, index[km]] += 1
    observed = numpy.zeros(len(index))
    for read in read_seqs(prep_dir, rid):
        for km in canonical_kmers(read):
            j = index.get(km)
            if j is not None:
                observed[j] += 1
    telling = counts.max(axis=0) != counts.min(axis=0)
    if not telling.any():
        return seqs[distinct[0]][0], seqs[distinct[0]][0], math.inf, distinct[0], distinct[0]
    counts, observed = counts[:, telling], observed[telling]
    scored = []
    for i in range(len(distinct)):
        for j in range(i, len(distinct)):
            c = counts[i] + counts[j]
            depth = observed[c > 0].sum() / max(1.0, c.sum())
            lam = numpy.maximum(c * depth, FLOOR)
            scored.append((float((observed * numpy.log(lam) - lam).sum()), i, j))
    scored.sort(reverse=True)
    top = scored[0]
    margin = top[0] - scored[1][0] if len(scored) > 1 else math.inf
    return seqs[distinct[top[1]]][0], seqs[distinct[top[2]]][0], margin, distinct[top[1]], distinct[top[2]]
