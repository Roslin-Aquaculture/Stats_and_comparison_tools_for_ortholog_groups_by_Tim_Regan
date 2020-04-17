#!/usr/bin/python
#Name/date: TimRegan/2020.04.15
#File: kofamscan_analyses.py
#Purpose: To analyse the highest scoring annotations from a Kofamscan output file
#Returns: A \t sep .csv with all KO definitions sorted by count, with cumulative score and avg. E-score for each hit. 
#Usgae:	python kofamscan_analyses.py <kofamscan_dir>

import sys
import glob
import pandas as pd

Kdir = sys.argv[1]
#The driectory given containing all Kofamscan output.txt files

path = r''+Kdir # use your path
all_files = glob.glob(path + "/*.fa.txt")
#Selects all .txt files in this directory - make sure all .txt are kofamscan outputs

def KofamscanStats(f):
	""" takes a .fa.txt Kofamscan output file and returns a .csv summary of the results """	
	df = pd.read_fwf(f, header=0, skiprows=[1], index_col=None)
	#This reads the whitespace delimited kofamscan output file, sets the headers and ignores the 2nd line of ---
	if 'KO definition' in df.columns:
		defs = df['KO definition'].value_counts().to_frame().reset_index()
		#Returns a df of KO definition Counts
		defs.columns = ['KO definition', 'Count']
		#Set Column titles
		defs = defs.sort_values('KO definition')
		#Sorted df (alphabetically) based on "KO definition". This matters for consistency when adding later columns. 
		Score = df.groupby('KO definition')['score'].sum()
		Eval = df.groupby('KO definition')['E-value'].mean()
		#This gives series of cumulative scores and avg. E-values for each KO def. 
		defs['Avg. E-value'] = Eval.values
		defs['Cum. Score'] = Score.values
		#This adds the series to the dataframe. 
		defs = defs.sort_values('Count', ascending=False)
		#This sorts the dataframe to have the highest count on top. 
		name = f.replace(".fa.txt", ".txt") 
		defs.to_csv(name, sep='\t', index=False)
	else:
		print "Error with "+f

def KofamscanThrshld(f):
	""" takes a .fa.txt Kofamscan output file and returns a .csv summary of the results WHICH MEET THRESHOLD """	
	df = pd.read_fwf(f, header=0, skiprows=[1], index_col=None)
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
			Score = is_thrshld.groupby('KO definition')['score'].sum()
			Eval = is_thrshld.groupby('KO definition')['E-value'].mean()
			#This gives series of cumulative scores and avg. E-values for each KO def. 
			defs['Avg. E-value'] = Eval.values
			defs['Cum. Score'] = Score.values
			#This adds the series to the dataframe. 
			defs = defs.sort_values('Count', ascending=False)
			#This sorts the dataframe to have the highest count on top. 
			name = f.replace(".fa.txt", ".Thrshld.txt") 
			defs.to_csv(name, sep='\t', index=False)
		else:
			return
	else:
		print "Error with Trshld for"+f

for filename in all_files:
	KofamscanStats(filename)
	KofamscanThrshld(filename)
