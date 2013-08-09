#!/usr/bin/python2.6
#Compare enhancer differences w/ motif location.
REFRESH=False #canonically false, but can be overwritten. THIS whole section should be parsed as args.
DIRECTORY='/net/home/carlesba/db/SSeq/enhancers/'
TARGETS= '/net/home/carlesba/project/dmel_targets'
COUNTS='/net/home/carlesba/db/SSeq/InsulinMotifs.count'
OUTFILE='/net/home/carlesba/db/SSeq/InsulinMotifs.tab'

import os.path
import sys
from math import exp
from math import log

# Must run script on its own / w/out "python"
if (len(sys.argv) > 1):
    DIRECTORY = sys.argv[1]
    TARGETS = sys.argv[2]
    COUNTS = sys.argv[3]
    OUTFILE = sys.argv[4]


def create_headers(filename,threshold):
    headers = []
    with open(filename,'r') as cfile:
        for line in cfile:
            line = line.split("\n")[0]
            line = line.split(" ")
            if (float(line[0]) > threshold):
                headers.append(line[1])
    return(set(headers))

HEADERS = create_headers(COUNTS,100) # NOTE That in fact, dm is only found at 30 instances. pnt at 100.

def gene_melt(gene,directory,headers,threshold):
    netin = directory + "/" + gene + ".net"
    motin = directory + "/" + gene + ".mot"
    if (os.path.exists(netin) and os.path.exists(motin)):
        # Find the length of the gene:
        length = sum(1 for line in open(netin))
        # Create dictionary of the file
        d={}
        with open(motin) as mot:
            for line in mot:
                line = line.split(" ")
                if float(line[2]) > threshold:
                    i = line[3] + line[4]
                    if not i in d:
                        d[i]=1 #also: if not i in d
                    else:
                        d[i]+=1
        # Create the out string:
        out = gene
        for name in headers:
            nameSY = name + "SY"
            nameNS = name + "NS"
            if not nameNS in d:
                NS = '0'
            else:
                NS = '%.3f' % (d[nameNS]/float(length))
            if not nameSY in d:
                SY = '0'
            else:
                SY = '%.3f' % (d[nameSY]/float(length))
            out = out + " " + SY + " " + NS
        out = out + " " + str(length) +"\n"
        return(out)
    else:
        return('EMPTY')
        # Only output dictionary names that exist already:
    
with open(TARGETS,'r') as tar:
    buffout = []
    for geneline in tar:
        g = geneline.split(" ")
        l = gene_melt(g[0],DIRECTORY,HEADERS,0.1)
        if not l == 'EMPTY':
            buffout.append(l)
with open(OUTFILE,'a') as out:
    # TODO Write headers first so that it looks like an R table.
    head = '"Gene" '
    for name in HEADERS:
        head = head + '"' + name + 'SY"' + ' "' + name + 'NS" '
    head = head + '"Length"\n'
    out.write(head)
    num =1
    for line in buffout:
        line = '"' + str(num) + '" ' + line
        out.write(line)
        num +=1
#Use occupancy values as da counts.



