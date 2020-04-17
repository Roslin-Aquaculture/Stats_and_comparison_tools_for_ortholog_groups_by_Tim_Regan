#!/usr/bin/python
#Name/date: TimRegan/2020.04.13
#File: Get_Bivalve_OGs.py
#Purpose: To ID and extract the names of orthogroups which have gene duplications retained in at least 50% bivalves analysed. 
#Usgae: python Get_Bivalve_OGs.py <Duplications.tsv> <Orthogroups.GeneCount.tsv> <Node>
from __future__ import division
#to enable decimal division
import sys
#to accept argv
import pandas as pd
import os
#for making new directories
from shutil import copyfile
#for copying .fasta files of Orthogroup sequences

dups = sys.argv[1]
#The Duplications.tsv file

OG_count = sys.argv[2]
#The counts per OG/species

node = str(sys.argv[3])
#This is the bivalve root node from the species tree. 

Bivalves = ["Bathymodiolus_platifrons","Crassostrea_gigas","Cristaria_plicata","Laternula_elliptica","Modiolus_philippinarum","Mya_arenaria","Mya_truncata","Mytilus_edulis","Mytilus_galloprovincialis","Pecten_maximus"]
Non_Bivalves =["Scutopus_ventrolineatus","Octopoteuthis_deletron","Octopus_bimaculoides","Vampyroteuthis_infernalis","Biomphalaria_glabrata","Lottia_gigantea","Lymnaea_stagnalis","Gymnomenia_pellucida","Wirenia_argentea","Laevipilina_hyalina","Acanthochitona_crinita","Gadila_tolmiei"]

def get_Node_OGs(f):
	#Load Duplications.tsv as dataframe, first row as headers
	df = pd.read_csv(f, sep="\t", index_col=None, header=0)
	##Get all OGs with <node> in second column, and support of 0.5 or greater in column 4. 
	is_node =  df['Species Tree Node']==node
	df_node = df[is_node]
	#Filtered file by bivalve node
	df_bv_support = df_node[df_node['Support'] >= 0.5]
	#Filtered by Support >= 0.5
	##Get avg. Support for every value of 'Orthogroup' and rank by that?
	##Instead rank by avg. genes/OG from Bv vs. avg. genes/OG from non-Bv 
	OGs = list(set(list(df_bv_support['Orthogroup'])))
	return OGs

def OG_FC(OG):
	avgbv = Bvdf.loc[OG,'mean']
	avgnb = NBdf.loc[OG,'mean']
	if avgnb == 0:
		FC = avgbv
	else:
		FC = avgbv/avgnb
	return FC

#Get orthogroups with duplications common to more than half bivalves ("OGs"). 
OGs = get_Node_OGs(dups)

#Get number of genes/OG and filter by Bivalves vs. Non Bivalves
countdf = pd.read_csv(OG_count, sep="\t", index_col="Orthogroup", header=0)
#Get average no. of genes/orthogroup for bivalves
Bvdf = countdf.loc[:, countdf.columns.isin(Bivalves)]
Bvdf['mean'] = Bvdf.mean(axis=1)
NBdf = countdf.loc[:, countdf.columns.isin(Non_Bivalves)]
NBdf['mean'] = NBdf.mean(axis=1)

#Make a corresponding list of corresponding fold change of Bv gene duplications for each OG 
#For bonus, copy .fa for each orthogroup to new directory
FC = []
os.mkdir('BivalveOGs')
for OG in OGs:
	FC.append(OG_FC(OG))
	fa = 'Results_Apr11/Orthogroup_Sequences/'+OG+'.fa'
	fac = 'BivalveOGs/'+OG+'.fa'
	copyfile(fa, fac)

#Make an output file with orthogroups and fold change of duplications in Bivalves vs. Non-Bivalves
outdf = pd.DataFrame(list(zip(OGs,FC)), columns=["Orthogroup","Fold Change Duplications"])

#Sort the dataframe based on fold change
foutdf = outdf.sort_values(by=["Fold Change Duplications"], ascending=False)

#Write this dataframe to a \t seperated .csv    
foutdf.to_csv('FC-Dups_in_Bivalve_Orthogroups.csv', sep="\t", index=False)
