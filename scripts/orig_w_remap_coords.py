#!/usr/bin/env python
"""Stream original VCF file with coordinates from remapped VCF output file.

NCBI remap has a bug where it dropsa secondary alleles in some cases:

1    6219287 G5566   TCACACA TCA,T   9797    PASS    LEN=4,6;TYPE=del,del    GT      1|2

1    6159227 G5566   TCACACA TCA     9797    PASS    LEN=4,6;TYPE=del,del;REMAP_ALIGN=FP     GT      1|2

This avoids this by writing the original record, but with the new coordinates.
"""
import sys

def main(orig_file, remap_file):
    remap_coords = read_remap_file(remap_file)
    with open(orig_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                parts = line.split("\t")
                if parts[2] in remap_coords:
                    new_chrom, new_start = remap_coords[parts[2]]
                    parts[0] = new_chrom
                    parts[1] = new_start
                    line = "\t".join(parts)
                else:
                    line = None
            if line:
                sys.stdout.write(line)

def read_remap_file(in_file):
    out = {}
    with open(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                chrom, start, vid = line.split("\t")[:3]
                out[vid] = (chrom, start)
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
