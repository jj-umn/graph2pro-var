#!/usr/bin/env python

from sys import argv

if len(argv) < 3:
	print('Usage: %s <input_file> <FDR value>' % argv[0])
	exit(0)

filename = argv[1]
FDRLevel = argv[2]
filename_FDR = argv[1]+"."+ FDRLevel +".tsv"
fw = open(filename_FDR,'w')
forward = 0.0
reverse = 0.0
titleHash = dict()
with open(filename) as file:
	for line in file:
		if line.startswith('#'):
			fw.write(line)
			continue
		lines = line.split('\t')
		title = lines[1]
		if title in titleHash :
			continue;
		else:
			titleHash[title] = 1
		evalue = float(lines[13])
		protein = lines[10]
		peptide = lines[9]
		peptide_new = peptide.replace('I','L')
		#if evalue <= 2.2E-11 :
		if protein.startswith('XXX_') :
			reverse += 1
		else:
			forward += 1
		spectrumQvalue = reverse * 1.0/forward
		fw.write(line)
		#print peptide, evalue
		if spectrumQvalue > float(FDRLevel) :
			break
fw.close()
print(forward,"\t",reverse,"\t",spectrumQvalue,"\t",evalue)

