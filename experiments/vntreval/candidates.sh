#!/bin/sh
# Candidate graphs of one aligner for every vntr-eval region package, for tools/evaluate_all.py:
#   sh experiments/vntreval/candidates.sh <method> <command> [flags...]
# runs `<command> [flags...] -r 1 <region>/hap32.fa` (plegma's or abPOA's MSA output) on every region of
# $VNTR_EVAL (~/plegma-bench/vntr-eval by default), turns each MSA into a graph with vntr-eval's msa_graph.py,
# the one every method is judged through, as candidates/<method>/<id>.gfa, and the seconds each took into
# work/plegma/<method>/seconds.txt. A region whose graph exists is skipped: name the method after the binary's
# commit (plegma_<commit>_<flags>) and each version's graphs and scores are made once.
method=$1
command=$2
shift 2
cd "${VNTR_EVAL:-$HOME/plegma-bench/vntr-eval}" || exit 1
mkdir -p candidates/$method work/plegma/$method
for d in regions/*/; do
  [ -s $d/hap32.fa ] || continue
  id=$(basename $d)
  [ -s candidates/$method/$id.gfa ] && continue
  start=$(date +%s.%N)
  timeout 3600 "$command" "$@" -r 1 $d/hap32.fa > work/plegma/$method/$id.msa.fa 2> work/plegma/$method/$id.err \
    || { echo "$id failed"; continue; }
  end=$(date +%s.%N)
  echo "$id $(echo "$end - $start" | bc)" >> work/plegma/$method/seconds.txt
  python3 -I tools/msa_graph.py work/plegma/$method/$id.msa.fa $d/hap32.fa candidates/$method/$id.gfa --quiet \
    > /dev/null || echo "$id graph failed"
done
