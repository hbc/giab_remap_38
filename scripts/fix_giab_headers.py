"""Fix GiaB headers

Add 'PLNCIIonWG' to header

and fixes alalts to:
    
##INFO=<ID=allalts,Number=.,Type=String
"""
import sys

wrote_extras = False
for line in sys.stdin:
    if line.find("ID=allalts") > 0:
        line = line.replace("Number=1,Type=Integer", "Number=.,Type=String")
    if line.startswith("##INFO") and not wrote_extras:
        wrote_extras = True
        new = '##INFO=<ID=PLNCIIonWG,Number=.,Type=String,Description="Genotype likelihoods for Ion">\n'
        sys.stdout.write(new)
    sys.stdout.write(line)

