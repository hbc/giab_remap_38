#!/bin/bash
# Prepare CrossMaped structural variant calls for structural variant calls
# Germline: Starting with a Genome in a Bottle deletion calls
set -eu -o pipefail
ref38=/cm/shared/apps/bcbio/20141204-devel/data/genomes/Hsapiens/hg38/seq/hg38.fa
python=python

# ## CrossMap
cm_dir=crossmap_grch38
chain=inputs/hg19ToHg38.over.chain
cm_out_base=GiaB_v2_19-38_crossmap

#
giab=inputs/giab-svclassify-deletions-2015-05-22.bed
syn4=inputs/synthetic_challenge_set4_tumour_25pctmasked_truth_sv_*.bed

outdir=sv_crossmap_hg38
mkdir -p $outdir
mkdir -p $outdir/ready

parallel  -k -j 1 -t "cat {} | sed 's/^\([0-9]\+\)\t/chr\1\t/g' | sed 's/^MT/chrM/g' | sed 's/^X/chrX/g' | sed 's/^Y/chrY/g' | CrossMap.py bed $chain /dev/stdin $outdir/{/.}-crossmap-hg38.bed" ::: $giab $syn4

parallel -k -j 1 -t "grep -v random {} | sort -V -k1,1 -k2,2n | bedtools merge -i - > $outdir/ready/{/}" ::: $outdir/*.bed
