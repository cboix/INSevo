#!/usr/bin/python2.6
REFRESH=False #canonically false, but can be overwritten. THIS whole section should be parsed as args.
DIRECTORY='/net/home/carlesba/db/SSeq/enhancers/'
TARGETS= '/net/home/carlesba/project/dmel_targets'

import os.path
import sys
from math import exp

# Must run script on its own / w/out "python"
if (len(sys.argv) > 1):
    DIRECTORY = sys.argv[1]
    TARGETS = sys.argv[2]

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
def read_PWMs(filename,S2expressed):
    with open(S2expressed,'r') as S2in:
        S2EX=[]
        names=[]
        for line in S2in:
            line = line.split("\n")[0]
            line = line.split(" ")
            S2EX.append(dict(name=line[0],fpkm=float(line[1])))
            names.append(line[0])
    with open(filename,'r') as infile:
        PWM=[]
        lines = []
        for line in infile:
            lines.append(line)
            if len(lines) >= 5:
                dic = read_motif_table(lines)
                if any(dic['fbgn'] in s for s in names):
                    PWM.append(dic)
                lines = []
        if len(lines) > 0:
            print('something left!')
        return(PWM)

def compare_motif(lines,pwm,consensus,threshold):
    # Could make this more efficient by cutting off a threshold earlier
    # To make this more efficient would be to have a matrix (length *4) which you can multiply against the pwm
    value=0
    breaker=0
    n=0
    header=lines[0].split(" ")
    last = float(header[2])
    for line in lines:
        line = line.split(" ")
        this = float(line[2])
        #if not continuous, break!
        if abs(this - last) > 1:
            breaker=1
            break
        #Compare the dmelanogaster sequence only.
        value = value + float(pwm[line[0].upper()][n])
        # This can and should be adjusted.
        if value < -6:
            breaker=1
            break
        last = this
        end = line[2]
        n=n+1
    o = 1/(1+exp(consensus - value))
    if (o >= threshold):
        d = dict(start=header[2],end=end,position=header[4],peaksep=header[3],v=value,occ=o,pwm=pwm['gene'],iden=pwm['iden'])
        if breaker==0: 
            return(d)
        else:
            return(0)
    else:
        return(0)

def find_motifs(filename,pwm,threshold):
    motif = []
    consensus = 0
    for i in range(0,len(pwm['T'])):
        consensus = consensus + max(float(pwm['A'][i]),float(pwm['G'][i]),float(pwm['C'][i]),float(pwm['T'][i]))
    with open(filename,'r') as enhfile:
        #with open('/home/carles/Labwork/StarrSeq/cbt.yak.div','r') as enhfile:
        lines = []
        for line in enhfile:
            lines.append(line)
            if len(lines) == pwm['length']:
                motif.append(compare_motif(lines,pwm,consensus,threshold))
                #Output score, and the start point -- so that I can plot them?
                #Append to list of dictionary and then select a threshold
                #   and here get rid of the earliest line:
                lines = lines[1:]
        # Select a threshold and return a certain amount of motifs.
        # ENDS here
    matches=[]
    for m in motif:
        if m != 0:
            matches.append(m)
    return(matches)

def compute_file(filename,PWM,threshold):
    motifs = []
    for pwm in PWM:
        out = find_motifs(filename,pwm,threshold)
        motifs.append(out)
    d = []
    for m in motifs:
        if len(m) > 0:
            for dic in m:
                d.append(dic)
    #processing to be done here to select motif threshold. 
    return(d)

def gene_compute(gene,directory,refresh):
    namein = directory + "/" + gene + '.yak.div'
    nameout = directory + "/" + gene + '.enh'
    if (os.path.exists(namein)):
        # If out file exists, don't do this unless refresh.
        if (not os.path.exists(nameout) or refresh):
            motifs=compute_file(namein,PWM,0.025)
            with open(nameout,'a') as out:
                for m in motifs:
                    outline = m['start'] + " " + m['end'] + " " + str(m['occ']) + " " + m['pwm'] + " " + m['iden'] + "\n"
                    out.write(outline)

PWM = read_PWMs('/net/home/carlesba/db/SSeq/PWMs','/net/home/carlesba/db/SSeq/S2mid')

# Loop over all gene targets and compute the enhancer searches for all:
with open(TARGETS,'r') as tar:
    for geneline in tar:
        g = geneline.split(" ")
        gene_compute(g[0],DIRECTORY,REFRESH)

#To test:
# pwm=PWM[1]
# motif=find_motifs('/home/carles/Labwork/StarrSeq/cbt.yak.div',pwm,0.25)
#start = time.time()
#motifs=compute_file('/home/carles/Labwork/StarrSeq/cbt.yak.div',PWM,0.25)
#elapsed = time.time() - start
