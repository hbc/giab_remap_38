## Genome in a Bottle NA12878 Human Genome 38 remapped validation set

Scripts to remap the [Genome in a Bottle](https://sites.stanford.edu/abms/giab) NA12878
validation variant calls to
[build 38 (GRCh38/hg38) of the human genome](http://genomeref.blogspot.co.uk/2013/12/announcing-grch38.html).

These convert the VCF calls and assessment region BED files from build 37 to build 38
coordinates using remapping. We take multiple remapping approaches for testing
purposes:

- [The NCBI remapping service](http://www.ncbi.nlm.nih.gov/genome/tools/remap)
- [CrossMap](https://crossmap.sourceforge.net) with [UCSC chain files](http://genome.ucsc.edu/cgi-bin/hgLiftOver)

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
