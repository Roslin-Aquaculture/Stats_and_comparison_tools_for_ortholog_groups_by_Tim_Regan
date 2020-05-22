#!/usr/bin/python
#Name/date: TimRegan/2020.04.29
#File: kofamscan_summaries.py
#version: 2.1
#Purpose: To take a dir. of Kofamscan analyses output files and return a single file summary of top hits
#Usgae:	python kofamscan_summaries.py <FC-Dups_in_Class_Orthogroups.txt> <kofamscan_dir> <number of top hits to display per OG>
#e.g. python kofamscan_summaries.py OGs_Results_Apr28 Gastropoda 5 

import sys
import glob
import pandas as pd
import os.path

Class = str(sys.argv[2]) 
OGf = str(sys.argv[1])+"/"+"FC-Dups_in_"+Class+"_Orthogroups.txt"
#This feeds the "FC-Dups_in_Class_Orthogroups.txt" OG summary file as a dataframe.
DupsDir = str(sys.argv[1])+"/"+Class+"_Dups"
#The driectory given containing all Kofamscan output.txt files (no other .txt files) 
TopHits = int(sys.argv[3])
#Select the number of top hits to summarise

def make_summaries(OG, Folder, Tophits):
	"""Takes an FC-OG.txt file, a folder of kofamscan summaries and the no. of KEGG hits to display per orthogroup"""
	finalcols = ["Orthogroup","Fold Change Duplications","Support","KO definition","Count","Avg. E-value","Avg. Score"]
	outdf = pd.DataFrame(columns = finalcols)
	outdfT = pd.DataFrame(columns = finalcols)
	"""This is the final output dataframe. Currently empty"""
	OGs = pd.read_csv(OG, sep='\t', index_col=None, header=0)
	stop = len(OGs.index)
	"""Read OG file as df and get number of orthogroups (stop)"""
	for i in range(stop):
		Og = OGs.loc[i].to_dict()
		outdf = outdf.append(Og, ignore_index=True, sort=False)
		outdfT = outdfT.append(Og, ignore_index=True, sort=False)
		"""This adds the first OG to a dictionary, then appends it to the final df""" 
		sumryAllf = Folder+"/All_hits/"+OGs.iloc[i,0] +".txt"
		sumryTrshldf = Folder+"/Threshold/"+OGs.iloc[i,0] +".Thrshld.txt"
		"""This is a path to the kofamscan summary file for this OG (orthogroup)"""
		if os.path.isfile(sumryAllf):
			orthodf = (pd.read_csv(sumryAllf, sep='\t', index_col=None, header=0)).loc[0:TopHits,:]
			outdf = outdf.append(orthodf, ignore_index=True, sort=False)
		else:
			pass	
		if os.path.isfile(sumryTrshldf):
			orthodfT = (pd.read_csv(sumryTrshldf, sep='\t', index_col=None, header=0)).loc[0:TopHits,:]
			outdfT = outdfT.append(orthodfT, ignore_index=True, sort=False)
		else:
			pass	
		"""If this summary file exists, add the first X lines of it to the final summary file. Else, next OG""" 
	outdf = outdf[finalcols]
	outdfT = outdfT[finalcols]
	fpath = Folder+'/'+Class+'_OGs0.5SuppDups_'+str(TopHits)+'KEGG.txt'
	outdf.to_csv(fpath, sep="\t", index=False)
	fTpath = Folder+'/'+Class+'_OGs0.5SuppDups_'+str(TopHits)+'KEGG_Thrshld.txt'
	outdfT.to_csv(fTpath, sep="\t", index=False)
	
make_summaries(OGf, DupsDir, TopHits)

#Sort descending filenames. strip .fa.txt
LostOGs = str(sys.argv[1])+"/"+Class+"_Lost_OGs.txt"
LostDir = str(sys.argv[1])+"/"+Class+"_Loss/"
SpecOGs = str(sys.argv[1])+"/"+Class+"_Specific_OGs.txt"
SpecDir = str(sys.argv[1])+"/"+Class+"_Specific/"

def make_summary_NonDup(OG, Folder, Tophits):
	"""Takes an FC-OG.txt file, a folder of kofamscan summaries and the no. of KEGG hits to display per orthogroup"""
	finalcols = ["Orthogroup","KO definition","Count","Avg. E-value","Avg. Score"]
	outdf = pd.DataFrame(columns = finalcols)
	outdfT = pd.DataFrame(columns = finalcols)
	"""This is the final output dataframe. Currently empty"""
	OGs = pd.read_csv(OG, sep='\t', index_col=None, header=None, usecols=[0], names=['Orthogroup'])
	stop = len(OGs.index)
	"""Read OG file as df and get number of orthogroups (stop)"""
	for i in range(stop):
		Og = OGs.loc[i].to_dict()
		outdf = outdf.append(Og, ignore_index=True, sort=False)
		outdfT = outdfT.append(Og, ignore_index=True, sort=False)
		"""This adds the first OG to a dictionary, then appends it to the final df""" 
		sumryAllf = Folder+"All_hits/"+OGs.iloc[i,0] +".txt"
		sumryTrshldf = Folder+"Threshold/"+OGs.iloc[i,0] +".Thrshld.txt"
		"""This is a path to the kofamscan summary file for this OG (orthogroup)"""
		if os.path.isfile(sumryAllf):
			orthodf = (pd.read_csv(sumryAllf, sep='\t', index_col=None, header=0)).loc[0:TopHits,:]
			outdf = outdf.append(orthodf, ignore_index=True, sort=False)
		else:
			pass	
		if os.path.isfile(sumryTrshldf):
			orthodfT = (pd.read_csv(sumryTrshldf, sep='\t', index_col=None, header=0)).loc[0:TopHits,:]
			outdfT = outdfT.append(orthodfT, ignore_index=True, sort=False)
		else:
			pass	
		"""If this summary file exists, add the first X lines of it to the final summary file. Else, next OG""" 
	outdf = outdf[finalcols]
	outdfT = outdfT[finalcols]
	sumtype = (Folder.split('_')[-1]).strip('/')
	fpath = Folder+Class+sumtype+'_OGs0.5Supp_'+str(TopHits)+'KEGG.txt'
	outdf.to_csv(fpath, sep="\t", index=False)
	fTpath = Folder+Class+sumtype+'_OGs0.5Supp_'+str(TopHits)+'KEGG_Thrshld.txt'
	outdfT.to_csv(fTpath, sep="\t", index=False)

make_summary_NonDup(SpecOGs, SpecDir, TopHits)
make_summary_NonDup(LostOGs, LostDir, TopHits)