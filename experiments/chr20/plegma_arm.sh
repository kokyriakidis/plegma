#!/bin/bash
# chr20's plegma arm and the genotyping test, from a pgrealign run's extract (the star arm's, so both arms realign the
# same regions):
#   bash experiments/chr20/plegma_arm.sh <work dir> <plegma binary> [plegma flags (-A)] [arm name (plegma)]
# <work dir>/star/chr20 holds pgrealign's prepare and extract (chr20.gbz, extract/) and its realigned graph; the
# plegma arm's MSAs (realign_plegma.py), spliced in by pgrealign replace, are checked by pgrealign verify; then the
# genotyping arms (test_arm.sh) for the original graph, pgrealign's and plegma's. A stage whose output exists is
# skipped.
set -u -o pipefail
W=$1
PLEGMA=$2
FLAGS=${3:--A}
ARM=${4:-plegma}
HERE=$(cd "$(dirname "$0")" && pwd)
PGREALIGN=${PGREALIGN:-$HOME/rivals-src/pgrealign}
KARAN=${KARAN:-$HOME/karan/build/cli/karan}
VG_SPLICE=${VG_SPLICE:-$HOME/.local/bin/vg}
export PATH=$HOME/.local/bin:$PATH
cd "$W" || exit 1
S=star/chr20
P=$ARM/chr20
mkdir -p $P
pg() { python3 -I -c "import sys; sys.path.insert(0, '$PGREALIGN'); from pgrealign.cli import main; sys.argv[0] = 'pgrealign'; sys.exit(main())" "$@"; }
step() { echo "[$(date +%T)] $*"; }

step "realign (plegma $FLAGS)"
if [ ! -s $P/msas.jsonl.gz ]; then
  /usr/bin/time -v python3 -I "$HERE/realign_plegma.py" "$PLEGMA" "$KARAN" $S/extract/packages.jsonl.gz $P/msas.jsonl.gz \
    --jobs 8 --timeout 600 --flags "$FLAGS" --resume > $P/realign.log 2>> $P/realign.err || { echo REALIGN_FAILED; exit 1; }
fi
grep -E "^done|Elapsed|Maximum resident" $P/realign.err
step replace
if [ ! -s $P/chr20.realigned.gbz ]; then
  pg replace --graph $S/chr20.gbz --extract $S/extract --msas $P/msas.jsonl.gz -o $P/chr20.realigned.gfa \
    --gbz $P/chr20.realigned.gbz --vg "$VG_SPLICE" -t 4 > $P/replace.log 2>&1 || { echo REPLACE_FAILED; tail $P/replace.log; exit 1; }
fi
step verify
pg verify --before $S/chr20.gbz --after $P/chr20.realigned.gbz -o $P/verify.json --vg "$VG_SPLICE" > $P/verify.log 2>&1
grep -E '"ok"|"paths_differ"' $P/verify.json
export VG=${VG:-$HOME/rivals-src/vg-out/vg-91d38c8} THREADS=${THREADS:-8}
# The original and star arms once (an earlier run's are kept), then this one.
for arm in original:$S/chr20.gbz star:$S/chr20.realigned.gbz $ARM:$P/chr20.realigned.gbz; do
  name=${arm%%:*}
  graph=${arm#*:}
  step "arm $name"
  if [ "$name" != "$ARM" ] && [ -s arm_$name/short.vcf.gz ]; then echo "arm $name: kept"; continue; fi
  bash "$HERE/test_arm.sh" "$W" arm_$name $graph > arm_$name.log 2>&1 || echo "arm $name failed"
  tail -4 arm_$name.log
done
touch ${ARM}_arm.done
