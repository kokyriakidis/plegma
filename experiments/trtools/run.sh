#!/bin/sh
# Every tandem-repeat tool on one vntr-eval region's distinct haplotypes, for summary.py:
#   run.sh <vntr-eval> <region id> <vamos catalogue .tsv.gz> <out dir>
# The tools come from $TRTOOLS/bin (a micromamba env with longdust, sdust, tantan, trf, ultra, vamos, minimap2 and
# samtools). vamos annotates aligned contigs: each haplotype is mapped (minimap2 asm20) to the region's CHM13 span
# and decomposed alone against the catalogue's motifs for the loci inside the span.
set -e
ve=$1
rid=$2
catalogue=$3
out=$4
B=${TRTOOLS:-$HOME/micromamba/envs/trtools}/bin
mkdir -p "$out"
python3 -I - "$ve/regions/$rid/hap32.fa" "$out/hap.fa" "$out/chm13.fa" <<'EOF'
import sys
names, seqs = [], []
for line in open(sys.argv[1]):
    line = line.strip()
    if line.startswith(">"):
        names.append(line[1:].split()[0]); seqs.append([])
    elif line:
        seqs[-1].append(line.upper())
seen = set()
with open(sys.argv[2], "w") as f, open(sys.argv[3], "w") as g:
    for n, s in zip(names, ("".join(x) for x in seqs)):
        if n.startswith("CHM13"):
            g.write(">chm13\n%s\n" % s)
        if s not in seen:
            seen.add(s)
            f.write(">%s\n%s\n" % (n.replace("#", "_"), s))
EOF
"$B/longdust" "$out/hap.fa" > "$out/longdust.bed"
"$B/longdust" -k8 -w20000 "$out/hap.fa" > "$out/longdust_k8.bed"
"$B/sdust" "$out/hap.fa" > "$out/sdust.bed"
"$B/tantan" -f4 -w 500 "$out/hap.fa" > "$out/tantan.tsv"
(cd "$out" && "$B/trf" hap.fa 2 7 7 80 10 50 2000 -h -ngs > trf.dat) || true
"$B/ultra" -p 500 -c --tsv --hs -o "$out/ultra.tsv" "$out/hap.fa" > /dev/null 2>&1
# vamos: the catalogue's loci inside the span, in the CHM13 span's coordinates.
python3 -I - "$ve/regions/$rid/region.json" "$catalogue" "$out/vamos.bed" <<'EOF'
import gzip, json, sys
r = json.load(open(sys.argv[1]))
lo, hi = r["span_start"], r["span_end"]
with gzip.open(sys.argv[2], "rt") as f, open(sys.argv[3], "w") as out:
    for line in f:
        c = line.split("\t")
        if c[0] == r["contig"] and int(c[2]) > lo and int(c[1]) < hi:
            out.write("chm13\t%d\t%d\t%s\n" % (max(0, int(c[1]) - lo + 1), min(hi - lo + 1, int(c[2]) - lo + 1), c[3]))
EOF
"$B/samtools" faidx "$out/chm13.fa"
"$B/samtools" faidx "$out/hap.fa"
mkdir -p "$out/vamos"
for name in $(grep ">" "$out/hap.fa" | cut -c2-); do
    "$B/samtools" faidx "$out/hap.fa" "$name" > "$out/vamos/$name.fa"
    "$B/minimap2" -ax asm20 --secondary=no "$out/chm13.fa" "$out/vamos/$name.fa" 2> /dev/null |
        "$B/samtools" sort -o "$out/vamos/$name.bam" - 2> /dev/null
    "$B/samtools" index "$out/vamos/$name.bam"
    "$B/vamos" --contig -b "$out/vamos/$name.bam" -r "$out/vamos.bed" -s "$name" -o "$out/vamos/$name.vcf" -t 1 \
        > /dev/null 2>&1 || true
    rm -f "$out/vamos/$name.fa"
done
