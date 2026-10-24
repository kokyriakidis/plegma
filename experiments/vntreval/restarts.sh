#!/bin/bash
# Restarts of plegma -A on vntr-eval regions: each with the first start's placement in a random order (-V's sixth
# number), to measure how far the search's local optima spread, against the default (seed 0).
#   bash experiments/vntreval/restarts.sh <vntr-eval> <plegma> <out dir> <seeds> <region> [<region> ...]
# Writes <out dir>/<region>.<seed>.msa.fa and <out dir>/summary.txt: per MSA, vntr-eval's unit cost over the
# optima (implied.py), plegma's objective (prices 2,2,1) and the affine sum of pairs.
set -u
VE=$1
PLEGMA=$2
OUT=$3
SEEDS=$4
shift 4
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
for r in "$@"; do
  for s in $(seq 0 $((SEEDS - 1))); do
    echo "$r $s"
  done
done | xargs -P "${JOBS:-2}" -n 2 sh -c "\"$PLEGMA\" -A -N -V -2,1,0,0,0,\$1 -r 1 \"$VE/regions/\$0/hap32.fa\" > \"$OUT/\$0.\$1.msa.fa\" 2> /dev/null"
for f in "$OUT"/*.msa.fa; do
  unit=$(python3 -I "$HERE/implied.py" "$VE" "$f" 2>&1 | grep -v "^\[" | sed "s/.*anchored\/optimal //; s/, induced.*//")
  more=$(python3 -I -c "import sys; sys.path.insert(0, '$HERE'); from diagnose import costs; u, r, a = costs('$f'); print('objective', u + r, 'affine', a)")
  echo "$(basename "$f") unit $unit $more"
done > "$OUT/summary.txt"
touch "$OUT/done"
