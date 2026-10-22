#!/usr/bin/env python3
"""Real haplotype pairs from the vntr-eval region packages, for build/seeding.

usage: python3 -I experiments/seeding/real_pairs.py <vntr-eval>/regions <out prefix> [random pairs per region] [seed]

For every region: CHM13 against every other distinct hap32 sequence, and up to `random pairs per region` pairs
of distinct non-CHM13 sequences. Writes <prefix>.fa (pair k: pk_a then pk_b) and <prefix>.meta.tsv: name,
region, stratum, period, copies, kind (chm13 or random), lengths. There is no true alignment: build/seeding
scores these against each pair's optimum.
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


def main():
    regions, prefix = sys.argv[1], sys.argv[2]
    per_region = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    k = 0
    with open(prefix + ".fa", "w") as fa, open(prefix + ".meta.tsv", "w") as meta:
        meta.write("name\tregion\tstratum\tperiod\tcopies\tkind\ta_len\tb_len\n")
        for region in sorted(os.listdir(regions)):
            path = os.path.join(regions, region)
            if not os.path.exists(os.path.join(path, "region.json")):
                continue
            info = json.load(open(os.path.join(path, "region.json")))
            records = read_fasta(os.path.join(path, "hap32.fa"))
            chm13 = [s for name, s in records if name.startswith("CHM13")]
            distinct = []
            for name, s in records:
                if s not in distinct:
                    distinct.append(s)
            others = [s for s in distinct if not chm13 or s != chm13[0]]
            pairs = [("chm13", chm13[0], s) for s in others] if chm13 else []
            randoms = [(a, b) for i, a in enumerate(others) for b in others[i + 1:]]
            rng.shuffle(randoms)
            pairs += [("random", a, b) for a, b in randoms[:per_region]]
            for kind, a, b in pairs:
                name = "p%d" % k
                fa.write(">%s_a\n%s\n>%s_b\n%s\n" % (name, a, name, b))
                meta.write("%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\n" % (name, region, info["stratum"], info.get("period"),
                                                              info.get("copies"), kind, len(a), len(b)))
                k += 1


if __name__ == "__main__":
    main()
