#!/bin/bash
# Prepare liftOver and remap 38 Genome in a Bottle reference
# starting with a standard GiaB 37 reference.
# Converts into chunks for processing then merged back into final
# VCF files.
set -eu -o pipefail
ref38=/cm/shared/apps/bcbio/20141204-devel/data/genomes/Hsapiens/hg38/seq/hg38.fa
python=python
version=v2_19

# ## CrossMap
cm_dir=crossmap_grch38
chain=inputs/hg19ToHg38.over.chain
cm_out_base=GiaB_v2_19-38_crossmap

mkdir -p $cm_dir
[[ -f $cm_dir/$cm_out_base-orig.vcf ]] || CrossMap.py vcf $chain inputs/GiaB_$version.vcf.gz $ref38 $cm_dir/$cm_out_base-orig.vcf
[[ -f $cm_out_base.vcf.gz ]] || cat $cm_dir/$cm_out_base-orig.vcf | $python scripts/map_GRCh38_hg38.py inputs/GRCh38_ensembl2UCSC.txt | $python scripts/fix_giab_headers.py | $python scripts/remove_duplicate_alleles.py | bcftools annotate -h inputs/GRCh38-contig-header.txt - | vt sort -w 10000000 -m full - | bgzip -c > $cm_out_base.vcf.gz
[[ -f $cm_out_base.vcf.gz.tbi ]] || tabix -f -p vcf $cm_out_base.vcf.gz
[[ -f $cm_dir/$cm_out_base-input.bed ]] || cat inputs/GiaB_${version}_regions.bed | sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g" > $cm_dir/$cm_out_base-input.bed

[[ -f $cm_dir/$cm_out_base-orig.bed ]] || CrossMap.py bed $chain $cm_dir/$cm_out_base-input.bed $cm_dir/$cm_out_base-orig.bed
[[ -f $cm_out_base-regions.bed ]] || cat $cm_dir/$cm_out_base-orig.bed | $python scripts/filter_bed_contigs.py $cm_out_base.vcf.gz | sort -k1,1 -k2,2n | bedtools merge -i - > $cm_out_base-regions.bed

# ## prepare split VCF and BED inputs
to_hg19_chroms='sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g"'
split_dir=source_split_grch37
mkdir -p $split_dir
[[ -f $split_dir/GIAB_001.vcf ]] || $python scripts/Split_Source_Create_Makefile.py inputs/GiaB_$version.vcf.gz
parallel "cat {} | vcf2bed.py | $to_hg19_chroms > {= s:.vcf:.bed: =}" ::: $split_dir/*.vcf

# ## remap

# variants
rm_dir=remap_grch38
mkdir -p $rm_dir
remap_outfile=GiaB_$version-38_remap.vcf.gz
parallel -k -j 1 -t "$python scripts/remap_remove_truncated.py {}  $rm_dir/{/.}.vcf" ::: $split_dir/*.vcf

parallel -k -j 5 -t "[[ -f $rm_dir/{/.}.vcf ]] || perl scripts/remap_api.pl --mode asm-asm --from GCF_000001405.13 --dest GCF_000001405.26 --annotation {}  --annot_out $rm_dir/{/.}.vcf --report_out $rm_dir/{/.}.37to38rpt.tsv" ::: $split_dir/*.vcf
parallel -k -j 1 --halt 2 -t "[[ -f {}.gz ]] || $python scripts/orig_w_remap_coords.py $split_dir/{/.}.vcf {} | $python scripts/map_GRCh38_hg38.py inputs/GRCh38_ensembl2UCSC.txt | $python scripts/fix_giab_headers.py | $python scripts/fix_nomatching_ref.py $ref38 | bcftools annotate -h inputs/GRCh38-contig-header.txt - | vt sort -w 10000000 -m full - | bgzip -c > {}.gz" ::: $rm_dir/*.vcf
parallel -k --halt 2 -t "[[ -f {}.tbi ]] || tabix -f -p vcf {}" ::: $rm_dir/*.vcf.gz
RM_INPUTS=`parallel -k 'echo "--variant {}"' ::: remap_grch38/GIAB_*.vcf.gz`
[[ -f $remap_outfile ]] || gatk-framework -Xms500m -Xmx4g -T CombineVariants -U ALL --suppressCommandLineHeader -R $ref38 --out $remap_outfile $RM_INPUTS

# bed file
rmr_dir=remap_grch38_regions
rmr_outfile=GiaB_$version-38_remap-regions.bed
mkdir -p $rmr_dir
split -l 100000 -d inputs/GiaB_${version}_regions.bed $rmr_dir/split_regions
parallel -k -j 7 -t "[[ -f {}-hg38.bed ]] || perl scripts/remap_api.pl --mode asm-asm --from GCF_000001405.13 --dest GCF_000001405.26 --annotation {}  --annot_out {}-hg38.bed --report_out {}-hg38.37to38rpt.tsv" ::: $rmr_dir/split_regions??

[[ -f $rmr_outfile ]] || cat $rmr_dir/split_regions*-hg38.bed | grep -v ^track | grep -v ^## | sort -k1,1 -k2,2n | bedtools merge -i - | $python scripts/map_GRCh38_hg38.py inputs/GRCh38_ensembl2UCSC.txt | $python scripts/filter_bed_contigs.py $remap_outfile > $rmr_outfile
[[ -f $rmr_dir/unmapped.bed ]] || $python scripts/find_nomap_from_remap.py $rmr_dir/unmapped.bed $rmr_dir/*-hg38.37to38rpt.tsv

# # liftover
lo_dir=liftover_grch38
chain=inputs/hg19ToHg38.over.chain
lor_outbase=$lo_dir/GiaB_$version-38_liftover
mkdir -p $lo_dir
[[ -f $lor_outbase.bed ]] || liftOver $cm_dir/$cm_out_base-input.bed $chain $lor_outbase.bed $lor_outbase.unmapped.txt
#parallel -k -j 1 -t "[[ -f $lo_dir/{/.}.bed ]] || liftOver {} $chain $lo_dir/{/.}.bed $lo_dir/{/.}.unmapped.txt" ::: $split_dir/*.bed

# ## Calculate non
$python scripts/create_nomap_37regions.py inputs/GiaB_${version}_regions.bed GiaB_${version}-37_prep_regions.bed $rmr_dir/unmapped.bed $lor_outbase.unmapped.txt

