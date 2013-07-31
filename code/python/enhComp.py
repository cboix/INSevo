#!/usr/bin/python2.6
#Compare enhancer differences w/ motif location.
REFRESH=False #canonically false, but can be overwritten. THIS whole section should be parsed as args.
DIRECTORY='/net/home/carlesba/db/SSeq/enhancers/'
TARGETS= '/net/home/carlesba/project/dmel_targets'

import os.path
import sys
from math import exp
from math import log

# Must run script on its own / w/out "python"
if (len(sys.argv) > 1):
    DIRECTORY = sys.argv[1]
    TARGETS = sys.argv[2]

# Same from before:
def read_motif_table(lines):
    first=lines[0].split("\n")[0]
    header=first.split("_")
    gene=header[0].split(">")
    fbgn=header[-1]
    d = dict(gene=gene[1],fbgn=fbgn,iden=first)
    for line in lines[1:5]:
        proc=line.split("\t")[:-1]
        name=proc[0].split(" ")[0]
        info=proc[1:]
        d.update(zip(name,[info]))
    d['length'] = len(d['T'])
    return(d)

#CREATION of the motif dictionaries:
def read_PWMs(filename):
    with open(filename,'r') as infile:
        PWM=[]
        lines = []
        for line in infile:
            lines.append(line)
            if len(lines) >= 5:
                PWM.append(read_motif_table(lines))
                lines = []
        if len(lines) > 0:
            print('something left!')
        return(PWM)

def motif_mutation(l,m,PWM):
    position = int(l[2]) - int(m['start'])
    #could be more efficient??
    pwm =[]
    for motif in PWM:
        if (motif['iden'] == m['iden']):
            pwm = motif
            break
    value = 0
    consensus = 0
    if (l[1] == '-' or l[1] == 'N'):
        o = 0
        #TODO CHANGE num fields
        out = " ".join(l) + " NS " + str(o) + " " + m['occ'] + " NA NA" + "\n"
    else:
        if (len(pwm) > 0):
            for i in range(0,len(pwm['T'])):
                shift = max(float(pwm['A'][i]),float(pwm['G'][i]),float(pwm['C'][i]),float(pwm['T'][i]))
                consensus = consensus + shift
                if (i != position):
                    value = value + shift
                else:
                    value = value + float(pwm[l[1].upper()][i])
            o = 1/(1 + exp(consensus - value)) 
            # Line to be printed:
            #HERE is a threshold that should be included in the results! A change of 2X in occupancy is reported.
            if (abs(log(o/float(m['occ']),2)) > 1):
                out = " ".join(l) + " NS " + str(o) + " " + m['occ'] + " " + m['gene'] + " " + m['iden'] + "\n"
            else: 
                out = " ".join(l) + " SY" + " NA NA NA NA" + "\n"
        else:
            out = " ".join(l) + " NS" + " NA NA NA NA" + "\n"
    return(out)

# Read file line by line, find differences, see if in motif or not. 
def compare_gene(gene,directory,PWM,refresh):
    namein = directory + "/" + gene + ".yak.div"
    motifin = directory + "/" + gene + ".enh"
    nameout = directory + "/" + gene + ".out"
    if (os.path.exists(namein) and os.path.exists(motifin)):
        if (not os.path.exists(nameout) or refresh):
            #Make motif matrix:
            motifs=[]
            buff=[]
            with open(motifin,'r') as motin:
                for line in motin:
                    line = line.split("\n")[0]
                    l = line.split(" ")
                    g = l[3]
                    iden = l[4]
                    motifs.append(dict(start=l[0],end=l[1],occ=l[2],gene=g,iden=iden))
            with open(namein) as f:
                #Go by line in f, instead of trying to cordon off the motifs first. 
                for line in f:
                    #Assign synonymous OR if in motif assign value:
                    line = line.split("\n")[0]
                    l = line.split(" ")
                    if (l[0].upper() != l[1].upper()):
                        # multiple in motif??
                        sy = 1;
                        for m in motifs:
                            if (int(m['end']) > int(m['start'])):
                                if (int(l[2]) <= int(m['end']) and int(l[2]) >= int(m['start'])):
                                    outline = motif_mutation(l,m,PWM)
                                    sy = 0
                                    break
                            else:
                                if (int(l[2]) >= int(m['end']) and int(l[2]) <= int(m['start'])):
                                    outline = motif_mutation(l,m,PWM)
                                    sy = 0
                                    break
                        if (sy == 1):
                            # Need to add extra fields to compensate for NS work:
                            outline = " ".join(l) + " SY" + " NA NA NA NA" + "\n"
                        buff.append(outline)
            with open(nameout,'a') as out:
                for w in buff: 
                    out.write(w)

PWM = read_PWMs('/net/home/carlesba/db/SSeq/PWMs')

with open(TARGETS,'r') as tar:
    for geneline in tar:
        g = geneline.split(" ")
        compare_gene(g[0],DIRECTORY,PWM,REFRESH)

# END Remove files? Other math?
