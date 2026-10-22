#!/bin/sh
# Theseus's candidate graphs for every vntr-eval region package, as candidates.sh makes plegma's:
#   sh experiments/vntreval/theseus.sh <theseus_msa> [method]
# theseus_msa -t 0 (its own scoring), at most 30 minutes a region. Theseus names its rows Sequence_1.. in input
# order: they are given back the input's names before msa_graph.py reads them.
binary=$1
method=${2:-theseus}
cd "${VNTR_EVAL:-$HOME/plegma-bench/vntr-eval}" || exit 1
mkdir -p candidates/$method work/$method
for d in regions/*/; do
  [ -s $d/hap32.fa ] || continue
  id=$(basename $d)
  [ -s candidates/$method/$id.gfa ] && continue
  start=$(date +%s.%N)
  timeout 1800 "$binary" -t 0 -s $d/hap32.fa -f work/$method/$id.raw.fa > work/$method/$id.log 2>&1 \
    || { echo "$id failed or timed out"; continue; }
  end=$(date +%s.%N)
  echo "$id $(echo "$end - $start" | bc)" >> work/$method/seconds.txt
  awk 'NR == FNR { if (/^>/) { sub(/^>/, ""); split($0, w, " "); names[++n] = w[1] } next }
       /^>/ { print ">" names[++k]; next } { print }' $d/hap32.fa work/$method/$id.raw.fa > work/$method/$id.msa.fa
  python3 -I tools/msa_graph.py work/$method/$id.msa.fa $d/hap32.fa candidates/$method/$id.gfa --quiet \
    > /dev/null || echo "$id graph failed"
done
