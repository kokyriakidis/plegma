#!/usr/bin/env python3
"""What each tandem-repeat tool says about one region's haplotypes, from run.sh's outputs.

usage: python3 -I experiments/trtools/summary.py <run.sh out dir> [<region.json>]

Per haplotype: its length; longdust's (k7 and k8 -w20000) longest low-complexity interval; sdust's intervals and
bases; tantan's, TRF's and ULTRA's longest repeat (its span, period, copies, and how many repeats it reports); vamos's
units and the distinct catalogue motifs it used. Then, per tool, the periods it reported across haplotypes.
"""
import collections
import json
import os
import sys


def fasta_lengths(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = 0
        elif name:
            out[name] += len(line)
    return out


def bed(path):
    out = collections.defaultdict(list)
    if os.path.exists(path):
        for line in open(path):
            c = line.split("\t")
            if len(c) >= 3:
                out[c[0]].append((int(c[1]), int(c[2])))
    return out


def longest(spans):
    return max(spans, key=lambda s: s[1] - s[0]) if spans else None


def main():
    d = sys.argv[1]
    region = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
    lengths = fasta_lengths(os.path.join(d, "hap.fa"))
    ld, ld8, sd = bed(os.path.join(d, "longdust.bed")), bed(os.path.join(d, "longdust_k8.bed")), bed(os.path.join(d, "sdust.bed"))
    tantan = collections.defaultdict(list)
    for line in open(os.path.join(d, "tantan.tsv")):
        c = line.split("\t")
        tantan[c[0]].append((int(c[1]), int(c[2]), int(c[3]), float(c[4])))
    trf = collections.defaultdict(list)
    name = None
    for line in open(os.path.join(d, "trf.dat")):
        if line.startswith("@"):
            name = line[1:].strip()
        elif name and line.strip():
            c = line.split()
            trf[name].append((int(c[0]) - 1, int(c[1]), int(c[2]), float(c[3]), c[13]))
    ultra = collections.defaultdict(list)
    for line in open(os.path.join(d, "ultra.tsv")):
        c = line.rstrip("\n").split("\t")
        if c[0] == "SeqID":
            continue
        ultra[c[0]].append((int(c[1]), int(c[2]), int(c[3]), float(c[6]), c[5]))
    vamos = {}
    for n in lengths:
        p = os.path.join(d, "vamos", n + ".vcf")
        if not os.path.exists(p):
            continue
        for line in open(p):
            if line.startswith("#"):
                continue
            info = dict(x.split("=", 1) for x in line.split("\t")[7].split(";") if "=" in x)
            units = info.get("ALTANNO_H1", "")
            labels = units.split("-") if units else []
            vamos[n] = (len(labels), len(set(labels)), len(info.get("RU", "").split(",")))
    if region:
        print("region %s: %s, annotated period %s, %s copies in CHM13, motif %s" % (
            region["region_id"], region["stratum"], region["period"], region["copies"], region["motif"][:60]))
    print("%-26s %6s | %-11s %-11s | %-9s | %-21s | %-21s | %-21s | %s" % (
        "haplotype", "length", "longdust", "ld k8", "sdust", "tantan span p copies", "TRF span p copies",
        "ULTRA span p copies", "vamos units/motifs"))

    def rep(rows):
        if not rows:
            return "-"
        r = max(rows, key=lambda x: x[1] - x[0])
        return "%d-%d %d %.0f (%d)" % (r[0], r[1], r[2], r[3], len(rows))

    periods = collections.defaultdict(collections.Counter)
    for n, L in lengths.items():
        a, b = longest(ld.get(n, [])), longest(ld8.get(n, []))
        s = sd.get(n, [])
        for tool, rows in (("tantan", tantan[n]), ("TRF", trf[n]), ("ULTRA", ultra[n])):
            if rows:
                periods[tool][max(rows, key=lambda x: x[1] - x[0])[2]] += 1
        v = vamos.get(n)
        print("%-26s %6d | %-11s %-11s | %-9s | %-21s | %-21s | %-21s | %s" % (
            n[:26], L, "%d-%d" % a if a else "-", "%d-%d" % b if b else "-",
            "%d/%dbp" % (len(s), sum(e - x for x, e in s)), rep(tantan[n]), rep(trf[n]), rep(ultra[n]),
            "%d/%d of %d" % v if v else "-"))
    for tool, c in periods.items():
        print("%s: longest repeat's period across haplotypes: %s" % (tool, dict(c.most_common())))
    cons = collections.Counter(max(r, key=lambda x: x[1] - x[0])[4] for r in ultra.values() if r)
    print("ULTRA consensus of the longest repeat: %s" % {k[:40]: v for k, v in cons.most_common(3)})
    cons = collections.Counter(max(r, key=lambda x: x[1] - x[0])[4] for r in trf.values() if r)
    print("TRF consensus of the longest repeat: %s" % {k[:40]: v for k, v in cons.most_common(3)})


if __name__ == "__main__":
    main()
