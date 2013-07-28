#!/usr/bin/python2.6
PREFIX='/net/home/carlesba' #Take from arguments in bash.
REFRESH=True #canonically false, but can be overwritten. THIS whole section should be parsed as args.
targets= PREFIX + '/project/dmel_targets'

from math import exp

def read_motif_table(lines):
    header=lines[0].split("_")[:-1]
    gene=header[0].split(">")
    d = dict(gene=gene[1])
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

def compare_motif(lines,pwm,consensus,threshold):
    #Could make this more efficient by cutting off a threshold earlier
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
        #Compare the dmelanogaster sequence only.
        value = value + float(pwm[line[0].upper()][n])
        last = this
        end = line[2]
        n=n+1
    o = 1/(1+exp(consensus - value))
    if (o >= threshold):
        d = dict(start=header[2],end=end,position=header[4],peaksep=header[3],v=value,occ=o)
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
                #output score, and the start point -- so that I can plot them?
                #append to list of dictionary and then select a threshold
                #and here get rid of the earliest line:
                lines = lines[1:]
        #select a threshold and return a certain amount of motifs.
        #ENDS here
    matches=[]
    for m in motif:
        if m != 0:
            matches.append(m)
    return(matches)


PWM = read_PWMs('/home/carles/Labwork/StarrSeq/PWMs')

#For testing purposes:
PWM2=PWM[:100]

def compute_file(filename):
    motifs = []
    for pwm in PWM2:
        out = find_motifs(filename,pwm,0.25)
        motifs.append(out)
    d = []
    for m in motifs:
        if len(m) > 0:
            for dic in m:
                d.append(dic)
    #processing to be done here to select motif threshold. 
    return(d)



#To test:
motif=find_motifs('/home/carles/Labwork/StarrSeq/cbt.yak.div',pwm,0.25)

motifs=compute_file('/home/carles/Labwork/StarrSeq/cbt.yak.div')


with open('/home/carles/Labwork/StarrSeq/cbt.enh','a') as out:
    for m in motifs:
        outline = m['start'] + " " + m['end'] + " " + m['position'] + " " + str(m['occ']) + " " + m['peaksep'] + "\n"
        out.write(outline)


dG=[]
for m in motif:
    if m != 0:
        dG.append(m['occ'])


occ=[]
for m in dG:
    occ.append(1/(1 + math.exp(consensus-m)))


        dG.append(1/(1+math.exp(m['v'])))

