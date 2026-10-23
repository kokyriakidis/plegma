#!/bin/bash
# plegma variants on a few vntr-eval regions, scored as vntr-eval scores them (implied.py), with time and memory:
#   bash experiments/vntreval/variants.sh <vntr-eval> <plegma> <out dir> <regions a,b,...> <label=flags> [...]
# Writes <out dir>/<label>.<region>.msa.fa and prints one line per variant and region: unit cost over the optima,
# seconds and peak memory. JOBS (2) runs that many at once.
set -u
VE=$1
PLEGMA=$2
OUT=$3
REGIONS=$4
shift 4
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
for v in "$@"; do
  for r in ${REGIONS//,/ }; do
    echo "$v|$r"
  done
done | xargs -P "${JOBS:-2}" -I{} bash -c '
  v="${1%%|*}"; r="${1##*|}"; label="${v%%=*}"; flags="${v#*=}"
  /usr/bin/time -f "%e %M" -o "'"$OUT"'/$label.$r.time" "'"$PLEGMA"'" $flags -r 1 "'"$VE"'/regions/$r/hap32.fa" > "'"$OUT"'/$label.$r.msa.fa" 2> /dev/null
' _ {}
for v in "$@"; do
  label="${v%%=*}"
  for r in ${REGIONS//,/ }; do
    unit=$(python3 -I "$HERE/implied.py" "$VE" "$OUT/$label.$r.msa.fa" 2>&1 | grep -v "^\[" | sed "s/.*anchored\/optimal //; s/, induced.*//")
    read secs mem < "$OUT/$label.$r.time"
    printf "%-14s %-9s unit %s  %7.1f s %6d MB\n" "$label" "$r" "$unit" "$secs" $((mem / 1024))
  done
done
