#!/usr/bin/env python
"""Fix calls that do not match the reference genome.
"""
import sys

from pyfaidx import Fasta


refdict = Fasta(sys.argv[1])

for line in sys.stdin:
    if line.startswith("#"):
        matches = True
    else:
        chrom, start, rid, refa = line.split("\t")[:4]
        start = int(start) - 1
        refbase = str(refdict[chrom][start:start + len(refa)])
        matches = (refbase == refa) and all(b.upper() in ["G", "A", "T", "C", "N"] for b in refbase)
    if matches:
        sys.stdout.write(line)
