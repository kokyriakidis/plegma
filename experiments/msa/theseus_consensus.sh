#!/bin/sh
# Theseus's consensus of one read set on stdout, as accuracy.py runs a command:
#   sh experiments/msa/theseus_consensus.sh <theseus_msa> <type> <reads.fa>
# type 2: the heaviest bundle (theseus); 3: weighted majority voting (theseus_mv), whose output's first record
# is the consensus, the column weights and gapped consensus after it. Default scoring, as the benchmark harness
# runs it.
binary=$1
kind=$2
reads=$3
out=$(mktemp)
"$binary" -t "$kind" -s "$reads" -f "$out" > /dev/null 2>&1
status=$?
awk 'NR == 1 || (NR > 1 && !/^>/ && !done) { print; next } /^>/ { done = 1 }' "$out" | head -2
rm -f "$out"
exit $status
