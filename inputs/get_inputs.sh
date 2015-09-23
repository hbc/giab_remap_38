#!/bin/bash
# Retrieve input GiaB and chain files for remapping
set -eu -o pipefail

# chain files

wget -c https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2UCSC.txt
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# small variants

wget -c -O GiaB_v2_19.vcf.gz ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/NISTIntegratedCalls_14datasets_131103_allcall_UGHapMerge_HetHomVarPASS_VQSRv2.19_2mindatasets_5minYesNoRatio_all_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs.vcf.gz
wget -c -O GiaB_v2_19_regions.bed.gz ftp://ftp.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/NISTv2.19/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.19_2mindatasets_5minYesNoRatio.bed.gz

tabix -f -p vcf GiaB_v2_19.vcf.gz
gunzip *.chain.gz
gunzip *.bed.gz

# structural variants
wget -O - ftp://ftp.ncbi.nih.gov/giab/ftp/technical/svclassify_Manuscript/Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed | grep -v ^Chr > giab-svclassify-deletions-2015-05-22.bed
wget -c https://s3.amazonaws.com/bcbio_nextgen/dream/synthetic_challenge_set4_tumour_25pctmasked_truth.tar.gz
tar -xzvpf synthetic_challenge_set4_tumour_25pctmasked_truth.tar.gz
