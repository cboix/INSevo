#!/usr/bin/python2.6
# Be sure to run w/ python2.6 (v2.6.5) instead of python (which is python v2.4.3)
# A sample of the code found in aacomp.py: 
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

#NOTE Touch the output files each time (rm in fact)
with open('/net/home/carlesba/project/divergence/Akt1.simdiv','r') as f:
    with open('/net/home/carlesba/project/divergence/Akt1.out','a') as out:
        for line in f:
            a = line.split(" ")
            a[-1] = a[-1][0]
            code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
            out.write(code)
