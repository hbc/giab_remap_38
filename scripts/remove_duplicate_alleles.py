#!/usr/bin/env python
"""Remove VCF lines with duplicated alleles that make GATK upset.

Also remove ambiguous alleles.
"""

import sys

for line in sys.stdin:
    passes = True
    if not line.startswith("#"):
        chr, pos, vid, ref, alts = line.split("\t", 7)[:5]
        alleles = [ref] + alts.split(",")
        passes = len(alleles) == len(set(alleles))
        if passes:
            letters = [x.upper() for x in ref + alts if x != ","]
            passes = len(set(letters) - set(["A", "G", "T", "C"])) == 0
    if passes:
        sys.stdout.write(line)
