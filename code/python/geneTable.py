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
    # Find the length of the gene:
    length = sum(1 for line in open(netin))
    # Create dictionary of the file
    d={}
    with open(motin) as mot:
        for line in mot:
            line = line.split(" ")
            if float(line[2]) > threshold:
                i = line[3] + line[4]
                if not d.has_key(i):
                    d[i]=1  #also: if not i in d
                else:
                    d[i]+=1
    # Create the out string:
    out = gene
    for name in headers:
        nameSY = name + "SY"
        nameNS = name + "NS"
        if not d.has_key(nameNS):
            NS = '0'
        else:
            NS = str(d[nameNS])
        if not d.has_key(nameSY):
            SY = '0'
        else:
            SY = str(d[nameSY])
        out = out + " " + SY + " " + NS
    # Only output dictionary names that exist already:

with open(TARGETS,'r') as tar:
    buffout = []
    for geneline in tar:
        g = geneline.split(" ")
        buffout.append(gene_melt(g[0],DIRECTORY,HEADERS,0.1))
with open(OUTFILE,'a') as out:
    # TODO Write headers first so that it looks like an R table.
    head = []
    for name in HEADERS:
        head = head + '"' + name + 'SY"' + ' "' + name + 'NS" '
    out.write(head)
    for line in buffout:
        out.write(line)

