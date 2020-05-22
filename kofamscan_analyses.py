#!/usr/bin/python
#Name/date: TimRegan/2020.04.15
#File: kofamscan_analyses.py
"""
Ensure kofamscan is downloaded. Make a python env for this. 
See: https://github.com/takaram/kofam_scan (or https://bioconda.github.io/recipes/kofamscan/README.html)
Ensure config.yml is updated (example here: https://taylorreiter.github.io/2019-05-11-kofamscan/)
#Follows on from kofamscan being run:

$load hmmer and ruby modules
module load anaconda
source activate /exports/eddie/scratch/tregan/Kofamscan/kofamscan
cd /exports/eddie/scratch/tregan/Kofamscan
ls </exports/eddie/scratch/tregan/Mollusc_peps/New_peps/OrthoFinder/OGs_Results_May16/Bivalvia_Specific>/*.fa > filestoprocess.txt
F=`sed -n ${SGE_TASK_ID}p < filestoprocess.txt`
kofamscan/bin/exec_annotation -o "$F".txt -f detail-tsv $F --tmp-dir=./tmp_"$f"
qsub -t 1-<No. of .fa files to run> Kofamscan_array.sh 
"""
#Purpose: To analyse the highest scoring annotations from a Kofamscan output file
#Returns: A \t sep .txt with all KO definitions sorted by count, with avg. score and E-score for each hit. 
#Usgae:	python kofamscan_analyses.py <kofamscan_dir> 
#or: for d in OGs_Results_May16/*/; do python kofamscan_analyses.py $d ; done

import sys
import os
import glob
import pandas as pd

Kdir = sys.argv[1]
#The driectory given containing all Kofamscan output.txt files

path = r''+Kdir # use your path
all_files = glob.glob(path + "/*.fa.txt")
#Selects all .txt files in this directory - make sure all .txt are kofamscan outputs

def KofamscanStats(f):
	""" takes a .fa.txt Kofamscan output file and returns a .csv summary of the results """	
	df = pd.read_csv(f, sep='\t', index_col=None, header=0, skiprows=[1])
	#This reads the whitespace delimited kofamscan output file, sets the headers and ignores the 2nd line of ---
	if 'KO definition' in df.columns:
		defs = df['KO definition'].value_counts().to_frame().reset_index()
		#Returns a df of KO definition Counts
		defs.columns = ['KO definition', 'Count']
		#Set Column titles
		defs = defs.sort_values('KO definition')
		#Sorted df (alphabetically) based on "KO definition". This matters for consistency when adding later columns. 
		Score = df.groupby('KO definition')['score'].mean()
		Eval = df.groupby('KO definition')['E-value'].mean()
		#This gives series of cumulative scores and avg. E-values for each KO def. 
		defs['Avg. E-value'] = Eval.values
		defs['Avg. Score'] = Score.values
		#This adds the series to the dataframe. 
		defs = defs.sort_values('Count', ascending=False)
		#This sorts the dataframe to have the lowest E-value on top. 
		name = f.replace(".fa.txt", ".txt") 
		name = name.replace(Kdir, (Kdir+"/All_hits/"))
		defs.to_csv(name, sep='\t', index=False)
	else:
		print "Error with "+f

def KofamscanThrshld(f):
	""" takes a .fa.txt Kofamscan output file and returns a .csv summary of the results WHICH MEET THRESHOLD """	
	df = pd.read_csv(f, sep='\t', index_col=None, header=0, skiprows=[1])
	#This reads the whitespace delimited kofamscan output file, sets the headers and ignores the 2nd line of ---
	if 'KO definition' in df.columns:
		is_thrshld = df[df['#'].notna()]
		#Filtered file by score>threshold as denoted by the first column by * (or simply, 'notna').
		if is_thrshld.empty == False:
			#Not all files will have hits which meet threshold. We want to ignore these
			defs = is_thrshld['KO definition'].value_counts().to_frame().reset_index()
			#Returns a df of KO definition Counts
			defs.columns = ['KO definition', 'Count']
			#Set Column titles
			defs = defs.sort_values('KO definition')
			#Sorted df (alphabetically) based on "KO definition". This matters for consistency when adding later columns. 
			Score = is_thrshld.groupby('KO definition')['score'].mean()
			Eval = is_thrshld.groupby('KO definition')['E-value'].mean()
			#This gives series of cumulative scores and avg. E-values for each KO def. 
			defs['Avg. E-value'] = Eval.values
			defs['Avg. Score'] = Score.values
			#This adds the series to the dataframe. 
			defs = defs.sort_values('Count', ascending=False)
			#This sorts the dataframe to have the lowest E-value on top. 
			name = f.replace(".fa.txt", ".Thrshld.txt") 
			name = name.replace(Kdir, (Kdir+"/Threshold/"))
			defs.to_csv(name, sep='\t', index=False)
		else:
			return
	else:
		print "Error with Trshld for"+f

ThreshDir = Kdir+'/'+"Threshold"
NonThreshDir = Kdir+'/'+"All_hits"
if not os.path.exists(ThreshDir):
	os.makedirs(ThreshDir)
if not os.path.exists(NonThreshDir):
	os.makedirs(NonThreshDir)

for filename in all_files:
	KofamscanStats(filename)
	KofamscanThrshld(filename)
