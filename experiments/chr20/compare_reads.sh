#!/bin/bash
# The same reads across chr20 arms: where their alignments differ (only reads near the realigned regions can) and
# which graph aligns them better.
#   bash experiments/chr20/compare_reads.sh <work dir> <arm> <arm> [<arm> ...]
# Per arm, every record of <work dir>/arm_<arm>/short.gaf.gz as (read: name and mate, from giraffe's fn/fp tags;
# mapping quality; matches; aligned length; alignment score), sorted by read into <work dir>/reads.<arm>.tsv; then
# experiments/chr20/compare_reads.py on the joined table.
set -eu -o pipefail
W=$1; shift
HERE=$(cd "$(dirname "$0")" && pwd)
export LC_ALL=C
for arm in "$@"; do
  out=$W/reads.$arm.tsv
  [ -s $out ] && continue
  gzip -dc $W/arm_$arm/short.gaf.gz | awk -F'\t' '{tag="s"; as=0
      for (i = 13; i <= NF; i++) { if ($i ~ /^fn:Z:/) tag = "1"; else if ($i ~ /^fp:Z:/) tag = "2"; else if ($i ~ /^AS:i:/) as = substr($i, 6) }
      if ($6 == "*") print $1 "/" tag "\t0\t0\t0\t0"; else print $1 "/" tag "\t" $12 "\t" $10 "\t" $11 "\t" as}' \
    | sort -S 6G --parallel=4 -k1,1 > $out.tmp && mv $out.tmp $out
done
first=$1; shift
cp $W/reads.$first.tsv $W/reads.joined.tsv
for arm in "$@"; do
  join -t $'\t' $W/reads.joined.tsv $W/reads.$arm.tsv > $W/reads.joined.tmp && mv $W/reads.joined.tmp $W/reads.joined.tsv
done
python3 -I $HERE/compare_reads.py $W/reads.joined.tsv $first "$@"
