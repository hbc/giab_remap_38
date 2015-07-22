#!/usr/bin/env python
"""Stream original VCF file with coordinates from remapped VCF output file.

NCBI remap has a bug where it dropsa secondary alleles in some cases:

1    6219287 G5566   TCACACA TCA,T   9797    PASS    LEN=4,6;TYPE=del,del    GT      1|2

1    6159227 G5566   TCACACA TCA     9797    PASS    LEN=4,6;TYPE=del,del;REMAP_ALIGN=FP     GT      1|2

This avoids this by writing the original record, but with the new coordinates.

For dual mapped variants, prefers ones that map to the standard chromosomes
used in the GiaB reference (1-22 + X).
"""
import sys

std_chroms = ["%s" % x for x in range(1, 23) + ["X"]]

def main(orig_file, remap_file):
    remap_coords = read_remap_file(remap_file)
    with open(orig_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                parts = line.split("\t")
                if parts[2] in remap_coords:
                    line = ""
                    for new_chrom, new_start in remap_coords[parts[2]]:
                        parts[0] = new_chrom
                        parts[1] = new_start
                        line += "\t".join(parts)
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
                if vid in out:
                    prev_chrom, prev_start = out[vid][0]
                    # new one is in standard and old was not, overwrite
                    if chrom in std_chroms and prev_chrom not in std_chroms:
                        out[vid] = [(chrom, start)]
                    # new one not in standard, ignore
                    elif prev_chrom in std_chroms and chrom not in std_chroms:
                        pass
                    # both not in standard, add the latest
                    elif prev_chrom not in std_chroms and chrom not in std_chroms:
                        out[vid].append((chrom, start))
                    # both in standard, add the latest
                    else:
                        out[vid].append((chrom, start))
                else:
                    out[vid] = [(chrom, start)]
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
