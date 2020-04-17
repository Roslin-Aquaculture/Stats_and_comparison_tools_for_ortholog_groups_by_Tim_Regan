#!/usr/bin/python
#Name/date: TimRegan/2020.04.16
#File: Summ_BvOG_KEGG.py
#Purpose: To provide a summary of the top KEGG hits and scores for each of the top Bivalve Orthogroups
#Use in directory with KEGG summaries from "kofamscan_analyses.py" output. 
#Usgae: python Summ_BvOG_KEGG.py < number of hits to display per Orthogroup>

import sys
#to accept argv
import pandas as pd
import os

hits = sys.argv[1]+1
#This is the number of hits per Orthogroup that will be displayed (ordered by count)

OGs = pd.read_csv("OrthogroupsN9Support0.8.txt", sep='\t', index_col=None, header=0)

sumcols = OGs.iloc[0,0] + ".txt"
sumcols = list((pd.read_csv(sumcols, sep='\t', header=0, index_col=None)).columns)
outcols = (list(OGs.columns)) + sumcols
#combine OGs.columns and sumcols
outdf = pd.DataFrame(columns = outcols)
#Make a new dataframe using columns from each file type

stop = len(OGs.index)
#Will run for all Orthogroups listed in OGs file. 
for i in range(stop):
	OG = OGs.loc[i].to_dict()
	#Takes first row of OGs and creates dict.
	outdf = outdf.append(OG, ignore_index=True)
	#Adds this dict. to outdf dataframe
	path = OGs.iloc[i,0] + ".txt"
	#From Orthogroup name, gets path for corresponding KEGG summary file
	orthodf = pd.read_csv(path, sep='\t', index_col=None, header=0)
	#Reads this KEGG summary file
	orthodf = orthodf.loc[0:hits,:]
	#Take first 5 hits
	outdf = outdf.append(orthodf, ignore_index=True)
outdf = outdf[outcols]
#Orders columns in a sensible fashion (they get muddled otherwise for some reason)

outdf.to_csv('OGsN9Supp0.8_5KEGG.txt', sep="\t", index=False)
#Writes a final output file. 
