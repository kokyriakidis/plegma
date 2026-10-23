#!/bin/sh
# HG002 chr20 short reads for the chr20 genotyping test: the public 30x NovaSeq PCR-free alignment to GRCh38
# (Google Brain's), chr20 fetched by range, as paired FASTQ (mates collated; a read whose mate lies off chr20 kept
# as a single). Production's own chr20 read set, which vntr-eval used, is not public.
#   sh experiments/chr20/reads.sh <out dir>
S=${SAMTOOLS:-$HOME/micromamba/envs/bench-phasers/bin/samtools}
U=https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x/HG002.novaseq.pcr-free.30x.dedup.grch38.bam
cd "$1" || exit 1
[ -s hg002_30x.bai ] || curl -s -o hg002_30x.bai $U.bai
$S view -b -F 0x900 -X $U hg002_30x.bai chr20 -o hg002.chr20.grch38.bam || exit 1
$S collate -u -O hg002.chr20.grch38.bam tmp_collate \
  | $S fastq -1 hg002.chr20.R1.fq.gz -2 hg002.chr20.R2.fq.gz -s hg002.chr20.single.fq.gz -0 /dev/null -n - || exit 1
ls -la hg002.chr20.*
