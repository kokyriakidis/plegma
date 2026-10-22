#!/usr/bin/env python3
"""Every tandem-repeat test region of vntr-eval and pgrealign, in one manifest.

usage: python3 -I experiments/regions/gather.py <vntr-eval repository> <pgrealign repository> <suite dir>

Writes `<suite dir>/manifest.tsv`, one line per region of each set:

| set | regions | from |
| :-- | :-- | :-- |
| `packages` | vntr-eval's region packages (regions/<id>/), five strata | regions/regions.tsv |
| `hardset_chr20` | vntr-eval's hard chr20 repeats, pair F1 under 0.8 after pgrealign | hardset/chr20.tsv |
| `vntr_ge1kb` | vntr-eval's benchmarked VNTRs of 1 kb or more | regions/vntr_ge1kb.bed |
| `catalog` | pgrealign's CHM13 v2.0 tandem-repeat catalogue, every class | catalogs/chm13v2.0_tr.bed.gz |

Columns: set, id, contig, start, end (CHM13 v2.0, 0-based half-open), class (the catalogue's STR, VNTR, SAT
or LC, else VNTR), period, motif, stratum (packages only), and `sequences`: the directory under `<suite dir>`
holding the region's haplotypes when they are on hand, else `-` (to be extracted from a graph, see
`extract.sh`). A package's directory links its `hap32.fa` (34 haplotypes: 32 HPRC, CHM13, GRCh38), `panel.fa.gz`
(every HPRC haplotype spanning it, from `hprc.fa.gz`) and, for evaluation only, `truth.fa` (HG002).
Also writes `<suite dir>/summary.tsv`: per set and class, how many regions, how many 1 kb or longer, and how many
have their sequences on hand.
"""
import collections
import csv
import gzip
import os
import sys


def packages(vntr_eval, suite):
    rows = []
    with open(os.path.join(vntr_eval, "regions", "regions.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            source = os.path.join(vntr_eval, "regions", r["region_id"])
            target = os.path.join(suite, "packages", r["region_id"])
            os.makedirs(target, exist_ok=True)
            for name, link in (("hap32.fa", "hap32.fa"), ("hprc.fa.gz", "panel.fa.gz"), ("truth.fa", "truth.fa")):
                path = os.path.join(source, name)
                there = os.path.join(target, link)
                if os.path.exists(path) and not os.path.lexists(there):
                    os.symlink(path, there)
            rows.append(("packages", r["region_id"], r["contig"], r["span_start"], r["span_end"], "VNTR",
                         r["period"], r["motif"], r["stratum"], os.path.join("packages", r["region_id"])))
    return rows


def hardset(vntr_eval):
    rows = []
    with open(os.path.join(vntr_eval, "hardset", "chr20.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            rows.append(("hardset_chr20", r["id"], r["contig"], r["ref_start"], r["ref_end"], r["classes"] or "VNTR",
                         r["period"], r["motif"], "-", "-"))
    return rows


def ge1kb(vntr_eval):
    rows = []
    with open(os.path.join(vntr_eval, "regions", "vntr_ge1kb.bed")) as f:
        for line in f:
            if line.startswith("#"):
                continue
            contig, start, end, rid, period, motif = line.rstrip("\n").split("\t")[:6]
            rows.append(("vntr_ge1kb", rid, contig, start, end, "VNTR", period, motif, "-", "-"))
    return rows


def catalog(pgrealign):
    rows = []
    with gzip.open(os.path.join(pgrealign, "catalogs", "chm13v2.0_tr.bed.gz"), "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            contig, start, end, rid, kind, period, motif = line.rstrip("\n").split("\t")[:7]
            rows.append(("catalog", rid, contig, start, end, kind, period, motif, "-", "-"))
    return rows


def main():
    vntr_eval, pgrealign, suite = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(suite, exist_ok=True)
    rows = packages(vntr_eval, suite) + hardset(vntr_eval) + ge1kb(vntr_eval) + catalog(pgrealign)
    with open(os.path.join(suite, "manifest.tsv"), "w") as f:
        f.write("set\tid\tcontig\tstart\tend\tclass\tperiod\tmotif\tstratum\tsequences\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    counts = collections.defaultdict(lambda: [0, 0, 0])
    for row in rows:
        key = (row[0], row[5])
        counts[key][0] += 1
        counts[key][1] += int(row[4]) - int(row[3]) >= 1000
        counts[key][2] += row[9] != "-"
    with open(os.path.join(suite, "summary.tsv"), "w") as f:
        f.write("set\tclass\tregions\tat_least_1kb\twith_sequences\n")
        for (name, kind), (n, long, have) in sorted(counts.items()):
            f.write("%s\t%s\t%d\t%d\t%d\n" % (name, kind, n, long, have))
    print(open(os.path.join(suite, "summary.tsv")).read(), end="")


if __name__ == "__main__":
    main()
