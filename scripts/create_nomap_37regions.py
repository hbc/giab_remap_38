#!/usr/bin/env python
"""Remove regions that do not liftover from the original 37 regions file.
"""
import collections
import sys

def main(orig_file, out_file, *nomap_files):
    nomap_regions = collections.defaultdict(int)
    for nomap_file in nomap_files:
        if nomap_file.find("remap") >= 0:
            nomap_regions = add_remap_file(nomap_file, nomap_regions)
        elif nomap_file.find("liftover") >= 0:
            nomap_regions = add_liftover_file(nomap_file, nomap_regions)
        else:
            raise ValueError("Unexpected file: %s" % nomap_file)
    with open(orig_file) as in_handle:
        with open(out_file, "w") as out_handle:
            for line in in_handle:
                chrom, start, end = line.strip().split()
                key = (chrom, int(start), int(end))
                if nomap_regions[key] < 2:
                    out_handle.write(line)

def add_remap_file(in_file, regions):
    """Leave contigs as is (NCBI style) and remove 1 from start coordinate.
    """
    with open(in_file) as in_handle:
        for line in in_handle:
            chrom, start, end = line.strip().split()
            key = (chrom, int(start) - 1, int(end))
            regions[key] += 1
    return regions

def add_liftover_file(in_file, regions):
    """Remove chr from contig and do not change start/end coordinates.
    """
    with open(in_file) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                chrom, start, end = line.strip().split()
                key = (chrom.replace("chr", ""), int(start), int(end))
                regions[key] += 1
    return regions

if __name__ == "__main__":
    main(*sys.argv[1:])
