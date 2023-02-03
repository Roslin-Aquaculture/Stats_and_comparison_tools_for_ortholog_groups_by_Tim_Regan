#!/usr/bin/python
#Name/date: TimRegan/2020.05.01
#File: Find_containing.py
#Version: 1
#Purpose: Finds files matching a given string, or list of strings.
#Usgae: python Find_containing.py <Directory of .fa files string> <file with gene names> 

"""e.g. C.gigas IAP family genes (grep IAP Crassostrea_gigas.fa > C.gigas_IAPs.txt) (ALSO BCL2)
48 IAPs annotated from the genome of Crassostrea,5 and further supported by the 61 found in Pinctada in this study. The S. glomerata genome encodes 80 IAP
Unique to bivalve.
Feed this file of gene names to this script as second argument.
"""

import sys
import os
import glob
from Bio import SeqIO

Fadir = sys.argv[1]
#The driectory given containing all fasta files
path = r''+Fadir # use your path
all_files = glob.glob(path + "/*.fa")

def match_pattern(patterns, fa):
	"""Requires dictionary with species specific header patterns in each sequence of fasta file
	returns a fasta file containing all sequneces specific to a given species from a given fasta"""
	count = 0
	for seq_record in SeqIO.parse(fa, "fasta"):
		for pattern in patterns:
			if pattern in seq_record.description:
				count += 1
	return count

#Takes as second argument of command a list of sequence names/headers i.e. gene names of interest. 
patterns = str(sys.argv[2])
pattern_list = []
with open(patterns, 'r') as patternf:
	for line in patternf:
		if line.startswith('>'):
			line = line.split()
			genename = line[0].strip('>')
			pattern_list.append(genename)

matchnames = "Match_counts"+patterns
with open(matchnames, 'a') as matchnames:
	for filename in all_files:
		counts = match_pattern(pattern_list, filename)
		if counts > 0:
			counts = str(counts)
			matchnames.write(filename+'\t'+counts+'\n')
		else:
			pass
