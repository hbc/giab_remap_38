"""Find input regions that don't convert over to new coordinates using remap.
"""
import sys

def main(out_file, *in_files):
    with open(out_file, "w") as out_handle:
        for in_file in in_files:
            with open(in_file) as in_handle:
                for line in in_handle:
                    parts = line.split("\t")
                    if not line.startswith("#") and parts[4] == "NULL":
                        chrom = parts[3]
                        start = parts[7]
                        end = parts[8]
                        out_handle.write("%s\t%s\t%s\n" % (chrom, start, end))

if __name__ == "__main__":
    main(*sys.argv[1:])
