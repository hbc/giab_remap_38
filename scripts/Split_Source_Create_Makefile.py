# coding: utf-8

# #Split_Source_Create_Makefile
# 
# take the big source VCF from GIAB, split into files with 10000 records and then create a makefile to run remap
# 

from __future__ import division
import sys
import os
import vcf
import csv
import math
import gzip
import datetime


source_vcf = sys.argv[1]
out_dir="source_split_grch37/"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

nist_id=1
file_ct=1
errs=0
header=[]
data_line=[]
try:
    with gzip.open(source_vcf, 'r') as infile:
        data=csv.reader(infile, delimiter="\t")
        for line in data:
            if line[0].startswith("#"):
                #save header to reprint each time
                header.append(line)
            else:
                #save the data lines
                data_line.append(line)
except IOError:
    print "Can't open %s" % source_vcf
    
num=len(data_line)/20000
fi_num=math.ceil(num)
ln_ct=0 
var_ct=0
out_file=""
bed_file=""
file_list=[]
for li in data_line:
    if ln_ct == 0:
        out_file="%sGIAB_%03d.vcf" % (out_dir, file_ct)
        file_list.append(out_file)
        out=open(out_file, 'w')
        for h in header:
            out.write("%s\n" % "\t".join(h))
    if ln_ct<20001:
        ln_ct += 1
        var_ct += 1
        li[2]="G%d" % var_ct
        out.write("%s\n" % "\t".join(li))
        
    else:
        out.close()
        ln_ct = 0
        file_ct += 1
        var_ct += 1
        out_file="%sGIAB_%03d.vcf" % (out_dir, file_ct)
        file_list.append(out_file)
        out=open(out_file, 'w')
        ln_ct += 1
        for h in header:
            out.write("%s\n" % "\t".join(h))
        li[2]="G%d" % var_ct
        out.write("%s\n" % "\t".join(li))
        
out.close()
