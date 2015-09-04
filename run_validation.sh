#!/bin/bash
#
# Evaluate GiaB comparison true and false positives with  vcfeval and hap.py
# for 37 and remapped/crossmapped 38.
#
set -eu -o pipefail

ref_dir=/cm/shared/apps/bcbio/20141204-devel/data/genomes/Hsapiens
bcbio_dir=bcbio_calls
out_dir=giab-hg38-validation-results

ref37=$ref_dir/GRCh37/seq/GRCh37.fa
ref37sdf=$ref_dir/GRCh37/rtg/GRCh37.sdf
ref38=$ref_dir/hg38/seq/hg38.fa
ref38sdf=$ref_dir/hg38/rtg/hg38.sdf

# 37

truth37=$bcbio_dir/input/GiaB_v2_19.vcf.gz
truth37r=$bcbio_dir/input/GiaB_v2_19-37_prep_regions.bed
comp37=$bcbio_dir/final/NA12878-3/NA12878-3-gatk-haplotype.vcf.gz

#hap.py --keep-scratch -R $truth37r -f $truth37r -r $ref37 -o GiaB_37 -V $truth37 $comp37 -l 1,2,3,4,5,6,7.8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X
[[ -f $out_dir/rtg/GiaB_37/done ]] || rtg vcfeval -b $truth37 --bed-regions $truth37r -c $comp37 -o $out_dir/rtg/GiaB_37 -t $ref37sdf

# 38 -- remap

r_truth38=$bcbio_dir/input/GiaB_v2_19-38_remap.vcf.gz
r_truth38r=$bcbio_dir/input/GiaB_v2_19-38_remap-regions.bed
comp38=$bcbio_dir/final/NA12878-1/NA12878-1-gatk-haplotype.vcf.gz
fb_comp38=$bcbio_dir/final/NA12878-1/NA12878-1-freebayes.vcf.gz

#hap.py -R $r_truth38r -f $r_truth38r -r $ref38 -o GiaB_38_remap -V $r_truth38 $comp38
[[ -f $out_dir/rtg/GiaB_38-remap/done ]] || rtg vcfeval -b $r_truth38 --bed-regions $r_truth38r -c $comp38 -o $out_dir/rtg/GiaB_38-remap -t $ref38sdf
[[ -f $out_dir/rtg/GiaB_38-remap-fb/done ]] || rtg vcfeval -b $r_truth38 --bed-regions $r_truth38r -c $fb_comp38 -o $out_dir/rtg/GiaB_38-remap-fb -t $ref38sdf

# 38 -- crossmap

c_truth38=$bcbio_dir/input/GiaB_v2_19-38_crossmap.vcf.gz
c_truth38r=$bcbio_dir/input/GiaB_v2_19-38_crossmap-regions.bed

[[ -f $out_dir/rtg/GiaB_38-crossmap/done ]] || rtg vcfeval -b $c_truth38 --bed-regions $c_truth38r -c $comp38 -o $out_dir/rtg/GiaB_38-crossmap -t $ref38sdf

# explore discordants
edir=check_discordants
mkdir -p $edir
hc_remap_fns=$edir/hc-remap-fns.vcf.gz
hc_remap_fns_37=$edir/hc-37-fns.vcf.gz
fb_remap_fns=$edir/fb-remap-fns.vcf.gz
[[ -f $hc_remap_fns ]] || bcftools isec $comp38 $out_dir/rtg/GiaB_38-remap/fn.vcf.gz -n =2 -w 1 -c all -O z -o $hc_remap_fns
[[ -f $hc_remap_fns_37 ]] || bcftools isec $comp37 $out_dir/rtg/GiaB_37/fn.vcf.gz -n =2 -w 1 -c all -O z -o $hc_remap_fns_37
[[ -f $fb_remap_fns ]] || bcftools isec $fb_comp38 $out_dir/rtg/GiaB_38-remap-fb/fn.vcf.gz -n =2 -w 1 -c all -O z -o $fb_remap_fns

# Debug issues with mapping quality
#[[ -f $edir/b37-mq.txt ]] || zcat $comp37 | bio-vcf --eval 'rec.info.mq' --set-header "mq" > $edir/b37-mq.txt
#[[ -f $edir/b38-mq.txt ]] || zcat $comp38 | bio-vcf --eval 'rec.info.mq' --set-header "mq" > $edir/b38-mq.txt
