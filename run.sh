#!/bin/bash
# Prepare liftOver and remap 38 Genome in a Bottle reference
# starting with a standard GiaB 37 reference.
# Converts into chunks for processing then merged back into final
# VCF files.
set -eu -o pipefail
ref38=/cm/shared/apps/bcbio/20141204-devel/data/genomes/Hsapiens/hg38/seq/hg38.fa
python=/cm/shared/apps/bcbio/20141204-devel/data/anaconda/bin/python
version=v2_19

# prepare split VCF and BED inputs
split_dir=source_split_grch37
mkdir -p $split_dir
[[ -f $split_dir/GIAB_001.vcf ]] || $python scripts/Split_Source_Create_Makefile.py inputs/GiaB_$version.vcf.gz
to_hg19_chroms='sed "s/^\([0-9]\+\)\t/chr\1\t/g" | sed "s/^MT/chrM/g" | sed "s/^X/chrX/g" | sed "s/^Y/chrY/g"'
parallel "cat {} | vcf2bed.py | $to_hg19_chroms > {= s:.vcf:.bed: =}" ::: $split_dir/*.vcf

# liftover
lo_dir=liftover_grch38
chain=inputs/hg19ToHg38.over.chain
mkdir -p $lo_dir
parallel -k -j 1 -t "[[ -f $lo_dir/{/.}.bed ]] || liftOver {} $chain $lo_dir/{/.}.bed $lo_dir/{/.}.unmapped.txt" ::: $split_dir/*.bed

# ## remap

# variants
rm_dir=remap_grch38
mkdir -p $rm_dir
remap_outfile=GiaB_$version-38_remap.vcf.gz
parallel -k -j 1 -t "[[ -f $rm_dir/{/.}.vcf ]] || perl scripts/remap_api.pl --mode asm-asm --from GCF_000001405.25 --dest GCF_000001405.26 --annotation {}  --annot_out $rm_dir/{/.}.vcf --report_out $rm_dir/{/.}.37to38rpt.tsv --gbench_out $rm_dir/{/.}.gbp" ::: $split_dir/*.vcf
parallel -k -j 1 --halt 2 -t "[[ -f {}.gz ]] || $python scripts/orig_w_remap_coords.py $split_dir/{/.}.vcf {} | $python scripts/map_GRCh38_hg38.py inputs/GRCh38_ensembl2UCSC.txt | $python scripts/fix_giab_headers.py | $python scripts/fix_nomatching_ref.py $ref38 | bcftools annotate -h inputs/GRCh38-contig-header.txt - | vt sort -w 10000000 -m full - | bgzip -c > {}.gz" ::: $rm_dir/*.vcf
parallel -k --halt 2 -t "[[ -f {}.tbi ]] || tabix -f -p vcf {}" ::: $rm_dir/*.vcf.gz
RM_INPUTS=`parallel -k 'echo "--variant {}"' ::: remap_grch38/GIAB_*.vcf.gz`
[[ -f $remap_outfile ]] || gatk-framework -Xms500m -Xmx4g -T CombineVariants -U ALL --suppressCommandLineHeader -R $ref38 --out $remap_outfile $RM_INPUTS

# bed file
rmr_dir=remap_grch38_regions
rmr_outfile=GiaB_$version-38_remap-regions.bed
mkdir -p $rmr_dir
split -l 100000 -d inputs/GiaB_${version}_regions.bed $rmr_dir/split_regions
parallel -k -j 1 -t "[[ -f {}-hg38.bed ]] || perl scripts/remap_api.pl --mode asm-asm --from GCF_000001405.25 --dest GCF_000001405.26 --annotation {}  --annot_out {}-hg38.bed --report_out {}-hg38.37to38rpt.tsv --gbench_out {}-hg38.gbp" ::: $rmr_dir/split_regions??

[[ -f $rmr_outfile ]] || cat $rmr_dir/split_regions*-hg38.bed | grep -v ^track | grep -v ^## | sort -k1,1 -k2,2n | bedtools merge -i - | $python scripts/map_GRCh38_hg38.py inputs/GRCh38_ensembl2UCSC.txt > $rmr_outfile
