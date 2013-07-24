#!/usr/bin/python2.6
PREFIX='/net/home/carlesba' #Take from arguments in bash.
CDorNC='cd' #another argument that can be taken in bash.
REFRESH=True #canonically false, but can be overwritten. THIS whole section should be parsed as args.
targets= PREFIX + '/project/dmel_targets'

import os.path
import sys
# Must run script on its own / w/out "python"
if (len(sys.argv) > 1):
    PREFIX = sys.argv[1]
    targets = sys.argv[2]
    CDorNC= sys.argv[3]

# Be sure to run w/ python2.6 (v2.6.5) instead of python (which is python v2.4.3)
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons,amino_acids))

def compacids(a,b):
    a = a.lower()
    b = b.lower()
    if ('-' in b):
        return "DEL"
    else:
        if codon_table[a] == codon_table[b]:
            return "SY"
        elif codon_table[b.lower()] == "*":
            return "NONSENSE"
        else:
            return "NS"

# Function to write the comparison for one gene
def genecomp(gene,prefix,cdnc,refresh):
    namef = prefix + '/project/divergence/' + gene + '.' + cdnc + '.simdiv'
    nameout = prefix + '/project/divergence/' + gene + '.' + cdnc + '.out'

	# If out file exists, don't do this unless refresh.
    if (not os.path.exists(nameout) or REFRESH):
		with open(namef,'r') as f:
		    with open(nameout,'a') as out:
		        for line in f:
		            a = line.split(" ")
		            a[-1] = a[-1][0]
		            code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
		            out.write(code)

# Loop over all gene targets
with open(targets,'r') as tar:
    for geneline in tar:
        g = geneline.split(" ")
        genecomp(g[0],PREFIX,CDorNC,REFRESH)
