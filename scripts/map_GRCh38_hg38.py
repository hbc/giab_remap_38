"""Map GRCh38 coordinates (1, 2, 3) to hg38 (chr1, chr2).

Takes one argument, the ensembl to UCSC mapping file. Reads from stdin, writes to stdout.

https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCh38_ensembl2UCSC.txt
"""
import sys

remap_file = sys.argv[1]
remap = {}
with open(remap_file) as in_handle:
    for line in in_handle:
        parts = line.strip().split()
        if len(parts) == 2:
            remap[parts[0]] = parts[1]

for line in sys.stdin:
    if not line.startswith("##contig"):
        if line.startswith("#"):
            sys.stdout.write(line)
        else:
            parts = line.split("\t")
            new_chrom = remap.get(parts[0])
            if new_chrom:
                parts[0] = new_chrom
                sys.stdout.write("\t".join(parts))
