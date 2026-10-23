#!/bin/bash
# The files vntr-eval's evaluate.py truth stage reads, laid out as in the vg-call-eval repository it points at
# (VNTR_EVAL_DIR), for every contig with a vntr-eval region:
#   bash experiments/vntreval/truth_setup.sh <eval dir> <vg-call-eval repository> <truth dir> <vntr-eval>
#   data/truth/CHM13v2.0_HG2-T2TQ100-V1.1_{stvar,smvar}.{vcf.gz,benchmark.bed}   GIAB v5.0q on CHM13 (T2T-Q100 v1.1)
#   work/wgs/<contig>/<contig>.fa                                                    CHM13 v2.0 (UCSC hs1)
#   work/wgs-mm095/<contig>/truth.<contig>.stvar.{vcf.gz,bed}                       the contig's part
#   work/wgs-mm095/score/<contig>.truth.norm.vcf.gz                                  split by vg-call-eval's own
#                                                                                    split_multiallelic, as bench_wgs
set -u -o pipefail
E=$1; VCE=$2; T=$3; R=$4
BCF=${BCFTOOLS:-$HOME/micromamba/envs/bench-phasers/bin/bcftools}
TABIX=${TABIX:-tabix}
P=CHM13v2.0_HG2-T2TQ100-V1.1
mkdir -p $E/data/truth $E/work/wgs $E/work/wgs-mm095/score
for kind in stvar smvar; do
  ln -sf $T/HG002_CHM13v2.0_v5.0q_$kind.vcf.gz $E/data/truth/${P}_$kind.vcf.gz
  ln -sf $T/HG002_CHM13v2.0_v5.0q_$kind.vcf.gz.tbi $E/data/truth/${P}_$kind.vcf.gz.tbi
  ln -sf $T/HG002_CHM13v2.0_v5.0q_$kind.benchmark.bed $E/data/truth/${P}_$kind.benchmark.bed
done
if [ ! -s $E/work/hs1.fa ]; then
  curl -s https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz | gzip -dc > $E/work/hs1.fa || exit 1
  samtools faidx $E/work/hs1.fa 2>/dev/null || $HOME/micromamba/envs/bench-phasers/bin/samtools faidx $E/work/hs1.fa
fi
SAM=$HOME/micromamba/envs/bench-phasers/bin/samtools
for c in $(cut -f3 $R/regions/regions.tsv | sed 1d | sort -u); do
  d=$E/work/wgs/$c; m=$E/work/wgs-mm095/$c; mkdir -p $d $m
  [ -s $d/$c.fa ] || { $SAM faidx $E/work/hs1.fa $c > $d/$c.fa && $SAM faidx $d/$c.fa; }
  if [ ! -s $m/truth.$c.stvar.vcf.gz ]; then
    $BCF view -r $c -Oz -o $m/truth.$c.stvar.vcf.gz $T/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz && $TABIX -f -p vcf $m/truth.$c.stvar.vcf.gz
    awk -v c=$c '$1 == c' $T/HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed > $m/truth.$c.stvar.bed
  fi
  out=$E/work/wgs-mm095/score/$c.truth.norm.vcf.gz
  [ -s $out ] || PATH=$(dirname $BCF):$PATH python3 -c "
import sys
sys.path.insert(0, '$VCE/scripts/tier2')
sys.path.insert(0, '$VCE/src')
from pathlib import Path
from truvari_sv import split_multiallelic
split_multiallelic(Path('$m/truth.$c.stvar.vcf.gz'), Path('$out'), Path('$d/$c.fa'))
" || echo "$c: normalisation failed"
  echo "$c ready"
done
