#!/usr/bin/python2.6
#Be sure to run w/ python2.6 (v2.6.5) instead of python (which is python v2.4.3)
bases = ["t", "c", "a", "g"]
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
codon_table = dict(zip(codons,amino_acids))

def compacids(a,b):
	a = a.lower()
	b = b.lower()
	if ("-" in b):
		return "DEL"
	else:
		if codon_table[a] == codon_table[b]:
			return "SY"
		elif codon_table[b.lower()] == "*":
			return "NONSENSE"
		else:
			return "NS"
with open("/net/home/carlesba/project/divergence/Akt1.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Akt1.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/CAP.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/CAP.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/chico.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/chico.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/dm.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/dm.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/eIF-4E.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/eIF-4E.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Flo-2.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Flo-2.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/foxo.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/foxo.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/gig.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/gig.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp1.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp1.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp2.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp2.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp3.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp3.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp4.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp4.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp5.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp5.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp6.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp6.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Ilp7.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Ilp7.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/InR.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/InR.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/lkb1.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/lkb1.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Pi3K21B.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Pi3K21B.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Pi3K92E.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Pi3K92E.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Pten.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Pten.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/raptor.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/raptor.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Rheb.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Rheb.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/rictor.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/rictor.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/RpS6.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/RpS6.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/S6k.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/S6k.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/sgg.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/sgg.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/sima.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/sima.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Sin1.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Sin1.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/SNF1A.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/SNF1A.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/SNF4Agamma.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/SNF4Agamma.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/step.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/step.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Thor.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Thor.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Tor.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Tor.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
with open("/net/home/carlesba/project/divergence/Tsc1.simdiv","r") as f:
	with open("/net/home/carlesba/project/divergence/Tsc1.out","a") as out:
		for line in f:
			a = line.split(" ")
			a[-1] = a[-1][0]
			code = " ".join(a) + " " + compacids(a[0],a[1]) + "\n"
			out.write(code)
