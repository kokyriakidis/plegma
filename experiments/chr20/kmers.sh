#!/bin/bash
# HG002's 29-mer counts from its chr20 short reads, as KFF for vg haplotypes: jellyfish, then vntr-eval's jf2kff
# (tools/c/jf2kff.cpp, checked there to sample identically to KMC's counts), as vntr-eval's full-graph test made them.
#   bash experiments/chr20/kmers.sh <work dir>
JF=${JELLYFISH:-$HOME/micromamba/envs/jellyfish/bin/jellyfish}
J2K=${JF2KFF:-$HOME/.local/bin/jf2kff}
cd "$1" || exit 1
$JF count -m 29 -C -s 500M -t 8 -o reads.jf <(gzip -dc hg002.chr20.R1.fq.gz hg002.chr20.R2.fq.gz hg002.chr20.single.fq.gz) || exit 1
$JF dump -c -L 2 reads.jf | $J2K HG002.chr20.kff 29 || exit 1
rm -f reads.jf
ls -la HG002.chr20.kff
