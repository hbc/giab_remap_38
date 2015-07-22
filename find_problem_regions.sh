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

out_dir=debug_fps
mkdir -p $out_dir

gunzip -c $remap_fps | bedtools coverage -counts -a $remap_regions -b - | awk '{if ($4 > 10 && ($3 - $2) > 500) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -r -k 5n > $out_dir/remap_fps.bed
gunzip -c $remap_fns | bedtools coverage -counts -a $remap_regions -b - | awk '{if ($4 > 10 && ($3 - $2) > 500) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -r -k 5n > $out_dir/remap_fns.bed
gunzip -c $cm_fps | bedtools coverage -counts -a $cm_regions -b - | awk '{if ($4 > 10 && ($3 - $2) > 500) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -r -k 5n > $out_dir/crossmap_fps.bed
gunzip -c $cm_fns | bedtools coverage -counts -a $cm_regions -b - | awk '{if ($4 > 10 && ($3 - $2) > 500) print $1,$2,$3,$4,$4 * 1000/($3 - $2)}' | sort -r -k 5n > $out_dir/crossmap_fns.bed
