#!/usr/bin/env python3
"""Read sets simulated from real VNTR haplotypes, for consensus accuracy in repeats.

usage: python3 -I experiments/msa/vntr_reads.py <vntr-eval>/regions <out dir> <regions per stratum> [seed]

For each chosen region (hotspot_vntr first, then the other strata), one hap32 haplotype is the truth; read sets
of N in {10, 30} full-length reads at uniform error rates in {0.01, 0.05, 0.10} (substitution, insertion and
deletion alike, no homopolymer bias) are written as <out>/vntr_<region>_N<N>_r<rate>.fa beside
<...>.truth.fa, in the layout experiments/msa/accuracy.py reads (a `sim_vntr/` directory).
"""
import json
import os
import random
import sys


def read_fasta(path):
    out = []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            out.append([line[1:].split()[0], ""])
        elif out:
            out[-1][1] += line.upper()
    return out


def reads(rng, truth, n, rate):
    out = []
    for _ in range(n):
        s = []
        for c in truth:
            r = rng.random()
            if r < rate / 3:
                s.append(rng.choice([x for x in "ACGT" if x != c]))
            elif r < 2 * rate / 3:
                continue
            elif r < rate:
                s.append(c)
                s.append(rng.choice("ACGT"))
            else:
                s.append(c)
        out.append("".join(s))
    return out


def main():
    regions, out_dir, per = sys.argv[1], sys.argv[2], int(sys.argv[3])
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    target = os.path.join(out_dir, "sim_vntr")
    os.makedirs(target, exist_ok=True)
    chosen = {}
    for region in sorted(os.listdir(regions)):
        path = os.path.join(regions, region, "region.json")
        if os.path.exists(path):
            stratum = json.load(open(path))["stratum"]
            chosen.setdefault(stratum, []).append(region)
    for stratum, names in sorted(chosen.items()):
        rng.shuffle(names)
        for region in names[:per]:
            haplotypes = [s for name, s in read_fasta(os.path.join(regions, region, "hap32.fa"))
                          if not name.startswith("CHM13") and not name.startswith("GRCh38")]
            truth = rng.choice(haplotypes)
            if len(truth) > 30000:
                continue
            for n in (10, 30):
                for rate in (0.01, 0.05, 0.10):
                    stem = os.path.join(target, "vntr_%s_%s_L%d_N%d_r%02d" % (stratum, region, len(truth), n, int(rate * 100)))
                    with open(stem + ".truth.fa", "w") as f:
                        f.write(">truth\n%s\n" % truth)
                    with open(stem + ".fa", "w") as f:
                        for k, s in enumerate(reads(rng, truth, n, rate)):
                            f.write(">r%d\n%s\n" % (k, s))


if __name__ == "__main__":
    main()
