## Genome in a Bottle NA12878 Human Genome 38 remapped validation set

Scripts to remap the [Genome in a Bottle](https://sites.stanford.edu/abms/giab) NA12878
validation variant calls to
[build 38 (GRCh38/hg38) of the human genome](http://genomeref.blogspot.co.uk/2013/12/announcing-grch38.html).

These convert the VCF calls and assessment region BED files from build 37 to build 38
coordinates using remapping. We take multiple remapping approaches for testing
purposes:

- [The NCBI remapping service](http://www.ncbi.nlm.nih.gov/genome/tools/remap)
- [CrossMap](https://crossmap.sourceforge.net) with [UCSC chain files](http://genome.ucsc.edu/cgi-bin/hgLiftOver)

## Results

- [Validation results](http://imgur.com/a/2Ezon)

- [Remapped truth sets and other files](http://biodata.s3-website-us-east-1.amazonaws.com/giab_hg38_remap/)
  - Genome in a Bottle regions for GRCh37 that map to build 38:
    GiaB_v2_19-37_prep_regions.bed

  - Crossmap hg38 liftover with UCSC chain files, regions and VCF file:
    GiaB_v2_19-38_crossmap-regions.bed,  GiaB_v2_19-38_crossmap.vcf.gz

  - NCBI remap hg38 regions and VCF file:
    GiaB_v2_19-38_remap-regions.bed, GiaB_v2_19-38_remap.vcf.gz

  - Validation results for 37 and two 38 methods using [rtg vcfeval](https://github.com/RealTimeGenomics/rtg-tools):
    giab-hg38-validation-results.tar.gz

### Usage

Download the inputs with:

    cd inputs && bash get_inputs.sh

Run the remapping with:

    bash run.sh

### Requirements

This depends on external tools to do the actual work:

- Python with [pyfaidx](https://github.com/mdshw5/pyfaidx)
- Perl with XML::XPath
- [CrossMap](https://crossmap.sourceforge.net) -- can be installed with conda
- GNU parallel
- bedtools
- bcftools
- [vt](https://github.com/atks/vt)
- [GATK MIT licensed scripts](https://github.com/chapmanb/gatk)
- vcf2bed.py from [vcflib](https://github.com/ekg/vcflib)  

The easiest way to install the Python dependencies is with 
[Miniconda](http://conda.pydata.org/miniconda.html). Then do:

    conda install -c bcbio crossmap pyfaidx

### Contributors

[Deanna Church](https://github.com/deannachurch)
[Brad Chapman](https://github.com/chapmanb/)
