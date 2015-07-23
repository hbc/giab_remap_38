#!/usr/bin/env python
"""Stream original VCF file with coordinates from remapped VCF output file.

NCBI remap has a bug where it dropsa secondary alleles in some cases:

1    6219287 G5566   TCACACA TCA,T   9797    PASS    LEN=4,6;TYPE=del,del    GT      1|2

1    6159227 G5566   TCACACA TCA     9797    PASS    LEN=4,6;TYPE=del,del;REMAP_ALIGN=FP     GT      1|2

This avoids this by writing the original record for multiple alleles, 
but with the new coordinates.

For dual mapped variants, prefers ones that map to the standard chromosomes
used in the GiaB reference (1-22 + X).
"""
import sys
import collections

std_chroms = ["%s" % x for x in range(1, 23) + ["X"]]

def main(orig_file, remap_file):
    remap_coords, header = read_remap_file(remap_file)
    with open(orig_file) as in_handle:
        for line in header:
            sys.stdout.write(line)
        for line in in_handle:
            if line.startswith("#"):
                line = None
            else:
                parts = line.split("\t")
                if parts[2] in remap_coords:
                    line = ""
                    for new_chrom, new_start, new_line in remap_coords[parts[2]]:
                        # if multiple secondary alleles, map the location
                        if len(parts[4].split(",")) > 1:
                            parts[0] = new_chrom
                            parts[1] = new_start
                            line += "\t".join(parts)
                        # otherwise, pass the remap line
                        else:
                            line += new_line
                else:
                    line = None
            if line:
                sys.stdout.write(line)

def read_remap_file(in_file):
    header = []
    out = collections.defaultdict(list)
    with open(in_file) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                header.append(line)
            else:
                chrom, start, vid = line.split("\t")[:3]
                if chrom in std_chroms:
                    out[vid].append((chrom, start, line))
    return out, header

if __name__ == "__main__":
    main(*sys.argv[1:])
