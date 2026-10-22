#!/usr/bin/env python3
"""pgrealign's alignment of one read set, as a consensus: its medoid star (star_msa, abPOA underneath), then each
column's most frequent letter, a gap winning ties (abPOA's -a 1, plegma's -a 1).

usage: python3 -I experiments/msa/pgrealign_consensus.py <pgrealign repository> <abpoa> <reads.fa>

Every read is an allele of weight 1 and there are no fragments; star_msa runs with realign's scores and
`max_cells`, its homology guard off (`--min-share 0`): at 10% error a read shares too few 15-mers with the centre
for the default 0.5, and the pipeline would leave the region as it was. Prints the consensus as FASTA.
"""
import collections
import sys


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            names.append(line[1:].split()[0])
            seqs.append([])
        elif line:
            seqs[-1].append(line)
    return names, ["".join(s) for s in seqs]


def main():
    repository, abpoa, path = sys.argv[1], sys.argv[2], sys.argv[3]
    sys.path.insert(0, repository)
    from pgrealign.realign import MAX_CELLS
    from pgrealign.star import Abpoa, star_msa

    names, seqs = read_fasta(path)
    alleles = [("r%d" % k, s, 1) for k, s in enumerate(seqs)]
    with Abpoa(exe=abpoa) as align:
        msa = star_msa(alleles, align=align, max_cells=MAX_CELLS)
    rows = [msa.rows[a[0]] for a in alleles]
    out = []
    for column in zip(*rows):
        counts = collections.Counter(column)
        best = "-"
        for letter, n in counts.items():
            if n > counts[best]:
                best = letter
        if best != "-":
            out.append(best)
    print(">Consensus_sequence\n%s" % "".join(out))


if __name__ == "__main__":
    main()
