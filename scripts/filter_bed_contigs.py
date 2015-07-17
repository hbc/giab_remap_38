#!/usr/bin/env python
"""Filter a BED file to only contain contigs present in the reference VCF file.
"""
import sys
import gzip

vcf_file = sys.argv[1]

vcf_contigs = set([])
with gzip.open(vcf_file) as in_handle:
    for line in in_handle:
        if not line.startswith("#"):
            vcf_contigs.add("%s\t" % line.split("\t", 2)[0])

vcf_contigs = tuple(vcf_contigs)
for line in sys.stdin:
    if line.startswith(vcf_contigs):
        sys.stdout.write(line)
