#!/usr/bin/env python3
"""Which pair of a graph's alleles a region's reads support best, by depth as well as presence: the information the
reads hold, whatever a genotyper makes of it.

Every canonical K-mer of the rows (flanked) is counted in each row and in the reads; a pair of rows (a, b) expects a
K-mer's read count at (count in a + count in b) times one haplotype's K-mer depth (the median read count of the
K-mers the flanks hold once, halved), plus a little for errors, and the pair whose Poisson likelihood of the observed
counts is highest is the reads' choice. A copy-number difference shows as depth, where read placement cannot see it.
"""
import math
import os

import numpy

K = 31
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
    # One haplotype's K-mer depth: the flanks' K-mers that every row holds once.
    flank = [index[km] for km in set(canonical_kmers(left)) | set(canonical_kmers(right))
             if km in index and all(counts[i, index[km]] == 1 for i in range(len(distinct)))]
    depth = max(1.0, float(numpy.median(observed[flank])) / 2) if flank else 10.0
    error = 0.02 * depth
    scored = []
    for i in range(len(distinct)):
        lam = (counts[i][None, :] + counts[i:]) * depth + error
        ll = (observed[None, :] * numpy.log(lam) - lam).sum(axis=1)
        for t, j in enumerate(range(i, len(distinct))):
            scored.append((float(ll[t]), i, j))
    scored.sort(reverse=True)
    top = scored[0]
    margin = top[0] - scored[1][0] if len(scored) > 1 else math.inf
    return seqs[distinct[top[1]]][0], seqs[distinct[top[2]]][0], margin, distinct[top[1]], distinct[top[2]]
