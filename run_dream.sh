#!/bin/bash
# Prepare CrossMapped hg38 VCF and region files for DREAM synthetic datasets
set -eu -o pipefail
ref38=/cm/shared/apps/bcbio/20141204-devel/data/genomes/Hsapiens/hg38/seq/hg38.fa
python=python

# ## CrossMap
cm_dir=crossmap_grch38
chain=inputs/hg19ToHg38.over.chain

syn4_base=synthetic_challenge_set4_tumour_25pctmasked_truth

mkdir -p $cm_dir
[[ -f $cm_dir/$syn4_base-orig.vcf ]] || CrossMap.py vcf $chain inputs/$syn4_base.vcf.gz $ref38 $cm_dir/$syn4_base-orig.vcf
[[ -f $syn4_base-crossmap-hg38.vcf.gz ]] || cat $cm_dir/$syn4_base-orig.vcf | $python scripts/map_GRCh38_hg38.py inputs/GRCh38_ensembl2UCSC.txt | $python scripts/fix_giab_headers.py | $python scripts/remove_duplicate_alleles.py | bcftools annotate -h inputs/GRCh38-contig-header.txt - | vt sort -w 10000000 -m full - | bgzip -c > $syn4_base-crossmap-hg38.vcf.gz
[[ -f $syn4_base-crossmap-hg38.vcf.gz.tbi ]] || tabix -f -p vcf $syn4_base-crossmap-hg38.vcf.gz
[[ -f $cm_dir/$syn4_base-input.bed ]] || cat inputs/${syn4_base}_regions.bed | sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g" > $cm_dir/$syn4_base-input.bed
[[ -f $syn4_base-regions-crossmap-hg38.bed ]] || cat $cm_dir/$syn4_base-input.bed | $python scripts/filter_bed_contigs.py $syn4_base-crossmap-hg38.vcf.gz | sort -V -k1,1 -k2,2n | bedtools merge -i - > $syn4_base-regions-crossmap-hg38.bed
