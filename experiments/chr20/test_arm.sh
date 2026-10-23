#!/bin/bash
# One arm of the chr20 genotyping test, as vntr-eval's work/full/chr20/test_arm_latest.sh runs one, short reads only:
#   bash experiments/chr20/test_arm.sh <work dir> <arm name> <chr20 graph .gbz>
# sample 32 HG002 haplotypes (vg haplotypes --include-reference, CHM13 and GRCh38 kept) by HG002's 29-mers from its
# chr20 short reads -> the CHM13 path must still spell chr20 -> giraffe indexes -> giraffe (paired, fragment
# 401.7 +- 166.4 as vntr-eval's stage 4, then the singles) -> GAF-Base -> vg call with production's flags. One vg
# binary throughout ($VG, vg 91d38c802 as vntr-eval's full-graph test). Needs <work dir>/HG002.chr20.kff (kmers.sh),
# <work dir>/hg002.chr20.{R1,R2,single}.fq.gz (reads.sh) and <work dir>/chr20.fa (CHM13 chr20).
set -u -o pipefail
W=$1; A=$2; G=$3
VG=${VG:-$HOME/rivals-src/vg-91d38c8/bin/vg}
GAFBASE=${GAFBASE:-$HOME/.local/bin/gaf-base}
THREADS=${THREADS:-4}
D=$W/$A; mkdir -p $D $W/tmp
export TMPDIR=$W/tmp
step() { echo "[$(date +%T)] $*"; }
TM() { /usr/bin/time -v "$@"; }
peak() { grep -E "Elapsed|Maximum resident" "$1" | tr -s " " | tr "\n" " "; echo; }

step "$A: sampling index"
[ -s $D/full.hapl ] || TM $VG autoindex -n -w sampling -G $G -p $D/full -t $THREADS -T $TMPDIR -M 16G > $D/sampling.log 2> $D/sampling.err || { echo "$A SAMPLING_INDEX_FAILED"; exit 1; }
peak $D/sampling.err
step "$A: sample 32 HG002 haplotypes"
[ -s $D/sampled.gbz ] || TM $VG haplotypes -v 2 -t $THREADS -i $D/full.hapl -k $W/HG002.chr20.kff --include-reference \
    --set-reference CHM13 --set-reference GRCh38 --num-haplotypes 32 -g $D/sampled.gbz $G > $D/haplotypes.log 2> $D/haplotypes.err || { echo "$A SAMPLING_FAILED"; exit 1; }
peak $D/haplotypes.err
$VG paths -M -x $D/sampled.gbz | awk 'NR>1{print $2, $3}' | sort | uniq -c
$VG paths -F -x $D/sampled.gbz -Q CHM13#0#chr20 | sed '1s/.*/>chr20/' > $D/chm13.fa
# chr20.fa is soft-masked: compare the letters, not their case.
if cmp -s <(grep -v '>' $D/chm13.fa | tr -d '\n' | tr a-z A-Z) <(grep -v '>' $W/chr20.fa | tr -d '\n' | tr a-z A-Z); then
  echo "$A: CHM13#0#chr20 identical to chr20.fa"
else
  echo "$A: CHM13 PATH DIFFERS"; exit 1
fi
rm -f $D/chm13.fa
step "$A: giraffe indexes"
[ -s $D/map.dist ] || TM $VG autoindex -n -w sr-giraffe -G $D/sampled.gbz -p $D/map -t $THREADS -T $TMPDIR -M 12G > $D/autoindex.log 2> $D/autoindex.err || { echo "$A AUTOINDEX_FAILED"; exit 1; }
peak $D/autoindex.err
step "$A: giraffe short reads"
if [ ! -s $D/short.gaf.gz ]; then
  TM $VG giraffe -Z $D/sampled.gbz -d $D/map.dist -m $D/map.shortread.withzip.min -z $D/map.shortread.zipcodes -o gaf -t $THREADS \
     -f $W/hg002.chr20.R1.fq.gz -f $W/hg002.chr20.R2.fq.gz --fragment-mean 401.7 --fragment-stdev 166.4 2> $D/giraffe.paired.err | bgzip -@ 2 > $D/short.paired.gaf.gz
  peak $D/giraffe.paired.err
  TM $VG giraffe -Z $D/sampled.gbz -d $D/map.dist -m $D/map.shortread.withzip.min -z $D/map.shortread.zipcodes -o gaf -t $THREADS \
     -f $W/hg002.chr20.single.fq.gz 2> $D/giraffe.single.err | bgzip -@ 2 > $D/short.single.gaf.gz
  (gzip -dc $D/short.paired.gaf.gz; gzip -dc $D/short.single.gaf.gz) | grep -v '^@' | bgzip -@ 4 > $D/short.gaf.gz && rm -f $D/short.paired.gaf.gz $D/short.single.gaf.gz
fi
gzip -dc $D/short.gaf.gz | awk -F'\t' '{n++; if($6=="*") u++} END{print "short alignments", n, "unmapped", u+0}'
step "$A: GAF-Base and call"
rm -f $D/fifo; mkfifo $D/fifo
$GAFBASE sort $D/short.gaf.gz -o $D/fifo -p > $D/short.sort.log 2>&1 &
TM $GAFBASE construct $D/fifo -r $D/sampled.gbz -o $D/short.gaf.db --overwrite > $D/short.gafdb.log 2> $D/short.gafdb.err || { echo "$A GAFDB_FAILED"; exit 1; }
wait; rm -f $D/fifo
TM $VG call $D/sampled.gbz -p CHM13#0#chr20 -d 2 -t $THREADS -s HG002 --progress --mismap-max 0.95 --read-likelihood --phased \
   --mosaic-out $D/short.mosaic.tsv --gaf-base $D/short.gaf.db --gbz-base $D/sampled.gbz > $D/short.vcf 2> $D/short.call.err || { echo "$A CALL_FAILED"; exit 1; }
peak $D/short.call.err
bgzip -f $D/short.vcf && tabix -f -p vcf $D/short.vcf.gz
rm -f $D/map.dist $D/map.*.min $D/map.*.zipcodes
step "$A done"
