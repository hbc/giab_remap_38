#!/bin/bash
#
# Find remapped genome segments with large false positive or negative counts

set -eu -o pipefail

remap_regions=GiaB_v2_19-38_remap-regions.bed
cm_regions=GiaB_v2_19-38_crossmap-regions.bed

eval_base=giab-hg38-validation-results/rtg/GiaB_38
cm_fps=$eval_base-crossmap/fp.vcf.gz
cm_fns=$eval_base-crossmap/fn.vcf.gz
remap_fps=$eval_base-remap/fp.vcf.gz
remap_fns=$eval_base-remap/fn.vcf.gz
orig_fns=giab-hg38-validation-results/rtg/GiaB_37/fn.vcf.gz

out_dir=debug_fps
mkdir -p $out_dir

unique_remap_fns=$out_dir/remap-unique-fns.vcf.gz
python scripts/find_unique_fns.py $remap_fns $orig_fns remap_grch38/*.tsv | bgzip -c > $unique_remap_fns

gunzip -c $unique_remap_fns | bedtools coverage -counts -a $remap_regions -b - | awk '{if ($4 > 3 && ($3 - $2) > 10) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -k 5nr > $out_dir/remap-unique-fns.bed

gunzip -c $remap_fps | bedtools coverage -counts -a $remap_regions -b - | awk '{if ($4 > 5 && ($3 - $2) > 50) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -k 5nr > $out_dir/remap_fps.bed
gunzip -c $remap_fns | bedtools coverage -counts -a $remap_regions -b - | awk '{if ($4 > 5 && ($3 - $2) > 50) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -k 5nr > $out_dir/remap_fns.bed
gunzip -c $cm_fps | bedtools coverage -counts -a $cm_regions -b - | awk '{if ($4 > 5 && ($3 - $2) > 50) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -k 5nr > $out_dir/crossmap_fps.bed
gunzip -c $cm_fns | bedtools coverage -counts -a $cm_regions -b - | awk '{if ($4 > 5 && ($3 - $2) > 50) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -k 5nr > $out_dir/crossmap_fns.bed

