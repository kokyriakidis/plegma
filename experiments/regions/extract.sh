#!/bin/sh
# Every target's haplotypes, from an HPRC graph, by pgrealign's prepare, regions and extract, a contig at a time:
#   sh experiments/regions/extract.sh <graph.gbz> <pgrealign repository> <suite dir> [contig ...]
# (all of chr1..chr22, chrX, chrY by default). Writes <suite dir>/extract/<contig>/ (prepare's contig graph,
# snarls and reference index, regions.tsv.gz, and extract's packages.jsonl.gz with each region's alleles and their
# weights); a contig whose packages exist is skipped. The contig graph is removed afterwards to save disk.
gbz=$1
repo=$2
suite=$3
shift 3
contigs=${*:-"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"}
vg=${VG:-$HOME/.local/bin/vg}
# regions and extract find vg on PATH.
PATH=$(dirname $vg):$PATH
export PATH
pgrealign() {
  python3 -I -c "import sys; sys.path.insert(0, '$repo'); from pgrealign.cli import main; sys.argv[0] = 'pgrealign'; sys.exit(main())" "$@"
}
for contig in $contigs; do
  out=$suite/extract/$contig
  [ -s $out/packages.jsonl.gz ] && continue
  mkdir -p $out
  echo "$contig: prepare $(date)"
  pgrealign prepare --gbz $gbz --contig $contig -o $out --vg $vg -t 4 > $out/prepare.log 2>&1 || { echo "$contig: prepare failed"; continue; }
  pgrealign regions --contig $contig --targets $suite/targets_ge1kb.tsv --snarls $out/$contig.snarls.pb \
    --ref $out/$contig.ref.tsv.gz -o $out/regions.tsv.gz > $out/regions.log 2>&1 || { echo "$contig: regions failed"; continue; }
  pgrealign extract --graph $out/$contig.gbz --regions $out/regions.tsv.gz -o $out --ref-sample CHM13 \
    > $out/extract.log 2>&1 || { echo "$contig: extract failed"; continue; }
  rm -f $out/$contig.gbz
  echo "$contig: done $(date)"
done
