#!/usr/bin/python
#Name/date: TimRegan/2020.04.23
#File: Get_Class_OGs.py
#Purpose: To ID and extract the names of orthogroups which have gene duplications retained in at least 50% of Class analysed. 
#Edit "Support threshold" for Duplications if 0.5 not desired.
#Also outputs Class specific OGs, OGs whcih are in every species apart form this Class.
#Also gives a heatmap of no. of genes in OGs with Duplications supported in the Class selected for all species 
#Usgae: python Get_Class_OGs.py <Results_folder> <Node> <Class>
#e.g. 
#python Get_Class_OGs.py Results_May16 N9 Cephalopoda
#python Get_Class_OGs.py Results_May16 N7 Bivalvia
#python Get_Class_OGs.py Results_May16 N8 Gastropoda

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
#to enable decimal division
import sys
#to accept argv
import pandas as pd
import os
#for making new directories
from shutil import copyfile
#for copying .fasta files of Orthogroup sequences

#Taxa
Classes = ['Bivalvia','Caudofoveata','Cephalapoda','Gastropoda','Monoplacophora','Polyplacophora','Scaphopoda','Solenogastres']
Bivalvia = ['Pinctada_fucata', 'Pecten_maximus', 'Mytilus_galloprovincialis', 'Mytilus_edulis', 'Mya_truncata', 'Mya_arenaria', 'Modiolus_philippinarum', 'Mizuhopecten_yessoensis', 'Laternula_elliptica', 'Cristaria_plicata', 'Crassostrea_virginica', 'Crassostrea_gigas', 'Bathymodiolus_platifrons']
Caudofoveata = ['Scutopus_ventrolineatus']
Cephalopoda = ['Vampyroteuthis_infernalis', 'Octopus_vulgaris', 'Octopus_sinensis', 'Octopus_bimaculoides', 'Octopoteuthis_deletron']
Gastropoda = ['Pomacea_canaliculata', 'Lymnaea_stagnalis', 'Lottia_gigantea', 'Elysia_chlorotica', 'Biomphalaria_glabrata', 'Aplysia_californica', 'Achatina_fulica']
Monoplacophora = ['Laevipilina_hyalina']
Polyplacophora = ['Acanthochitona_crinita']
Scaphopoda = ['Gadila_tolmiei']
Solenogastres = ['Wirenia_argentea', 'Gymnomenia_pellucida']
Mollusca = Caudofoveata+Bivalvia+Cephalopoda+Gastropoda+Monoplacophora+Polyplacophora+Scaphopoda+Solenogastres

#Dictionary with each species and corresponding Class. 
SpeClaD = {'Acanthochitona_crinita': 'Polyplacophora', 'Achatina_fulica': 'Gastropoda', 'Aplysia_californica': 'Gastropoda', 'Bathymodiolus_platifrons': 'Bivalvia', 'Biomphalaria_glabrata': 'Gastropoda', 'Crassostrea_gigas': 'Bivalvia', 'Crassostrea_virginica': 'Bivalvia', 'Cristaria_plicata': 'Bivalvia', 'Elysia_chlorotica': 'Gastropoda', 'Gadila_tolmiei': 'Scaphopoda', 'Gymnomenia_pellucida': 'Solenogastres', 'Laevipilina_hyalina': 'Monoplacophora', 'Laternula_elliptica': 'Bivalvia', 'Lottia_gigantea': 'Gastropoda', 'Lymnaea_stagnalis': 'Gastropoda', 'Mizuhopecten_yessoensis': 'Bivalvia', 'Modiolus_philippinarum': 'Bivalvia', 'Mya_arenaria': 'Bivalvia', 'Mya_truncata': 'Bivalvia', 'Mytilus_edulis': 'Bivalvia', 'Mytilus_galloprovincialis': 'Bivalvia', 'Octopoteuthis_deletron': 'Cephalapoda', 'Octopus_bimaculoides': 'Cephalapoda', 'Octopus_sinensis': 'Cephalapoda', 'Octopus_vulgaris': 'Cephalapoda', 'Pecten_maximus': 'Bivalvia', 'Pinctada_fucata': 'Bivalvia', 'Pomacea_canaliculata': 'Gastropoda', 'Scutopus_ventrolineatus': 'Caudofoveata', 'Vampyroteuthis_infernalis': 'Cephalapoda', 'Wirenia_argentea': 'Solenogastres'}

node = str(sys.argv[2])
#This is the  root node from the species tree. 
Support_threshold = 0.5

TaxClass = str(sys.argv[3])
TaxList = eval(TaxClass)
NonClass = [s for s in Mollusca if s not in TaxList]

dups = str(sys.argv[1])+'/Gene_Duplication_Events/Duplications.tsv'
#The Duplications.tsv file
#Load in Duplications.tsv filtered by node of interest
dupf = pd.read_csv(dups, sep="\t", index_col=None, header=0)
is_node =  dupf['Species Tree Node']==node
dupf_node = dupf[is_node]

OG_count = str(sys.argv[1])+'/Orthogroups/Orthogroups.GeneCount.tsv'
#Load in this file (counts per OG/species)
countdf = pd.read_csv(OG_count, sep="\t", index_col="Orthogroup", header=0)
#Get average no. of genes/orthogroup for Class
Cldf = countdf.loc[:, countdf.columns.isin(TaxList)]
ClSpecDf = Cldf[(Cldf != 0).all(1)]
ClSpecDf = list(ClSpecDf.index.values)
ClSpecLoss = Cldf[(Cldf == 0).all(1)]
ClSpecLoss = list(ClSpecLoss.index.values)

NCdf = countdf.loc[:, countdf.columns.isin(NonClass)]
ClSpecNDf = ClSpecNDf = NCdf[(NCdf == 0).all(1)]
ClSpecNDf = list(ClSpecNDf.index.values)
NClSpecG = NCdf[(NCdf != 0).all(1)]
NClSpecG = list(NClSpecG.index.values)

ClassSpecOGs = set(ClSpecDf).intersection(ClSpecNDf)
ClassSpecOGs = list(sorted(ClassSpecOGs))
ClassSpecLoss = set(ClSpecLoss).intersection(NClSpecG)
ClassSpecLoss = list(sorted(ClassSpecLoss))

Cldf['mean'] = Cldf.mean(axis=1)
NCdf['mean'] = NCdf.mean(axis=1)


def get_Node_OGs(n):
	dupf_bv_support = dupf_node[dupf_node['Support'] >= n]
	#Filtered by Support >= n
	OGs = list(set(list(dupf_bv_support['Orthogroup'])))
	return OGs

def OG_Supp(OG):
	is_OG = dupf_node['Orthogroup']==OG
	dupf_OG = dupf_node[is_OG]
	Supp = dupf_OG['Support'].max()
	Supp = "{:.2f}".format(Supp)
	return Supp

def OG_FC(OG):
	avgbv = Cldf.loc[OG,'mean']
	avgnb = NCdf.loc[OG,'mean']
	if avgnb == 0:
		FC = float("inf")
	else:
		FC = avgbv/avgnb
		FC = "{:.2f}".format(FC)
	return FC

#Get orthogroups with duplications common to more than half of Class ("OGs"). 
OGs = get_Node_OGs(Support_threshold)

#Make a corresponding list of corresponding fold change and max Support value of Bv gene duplications for each OG 
#For bonus, copy .fa for each orthogroup to new directory
OGsDir = 'OGs_'+str(sys.argv[1])
if not os.path.exists(OGsDir):
    os.makedirs(OGsDir)
ClassOGsDupsDir = OGsDir + '/' + TaxClass + "_Dups"
if not os.path.exists(ClassOGsDupsDir):
    os.makedirs(ClassOGsDupsDir)
FC = []
Supp = []
for OG in OGs:
	FC.append(OG_FC(OG))
	Supp.append(OG_Supp(OG))
	fa = str(sys.argv[1])+'/Orthogroup_Sequences/'+OG+'.fa'
	fac = ClassOGsDupsDir+'/'+OG+'.fa'
	copyfile(fa, fac)

ClassOGsLossDir = OGsDir + '/' + TaxClass + "_Loss"
if not os.path.exists(ClassOGsLossDir):
    os.makedirs(ClassOGsLossDir)
for O in ClassSpecLoss:
	fa = str(sys.argv[1])+'/Orthogroup_Sequences/'+O+'.fa'
	fac = ClassOGsLossDir+'/'+O+'.fa'
	copyfile(fa, fac)

ClassOGsSpecDir = OGsDir + '/' + TaxClass + "_Specific"
if not os.path.exists(ClassOGsSpecDir):
    os.makedirs(ClassOGsSpecDir)
for O in ClassSpecOGs:
	fa = str(sys.argv[1])+'/Orthogroup_Sequences/'+O+'.fa'
	fac = ClassOGsSpecDir+'/'+O+'.fa'
	copyfile(fa, fac)

#Make an output file with orthogroups and fold change of duplications in TaxClass vs. Non-TaxClass
outdf = pd.DataFrame(list(zip(OGs,FC,Supp)), columns=["Orthogroup","Fold Change Duplications","Support"])
outdf['Fold Change Duplications'] = pd.to_numeric(outdf['Fold Change Duplications'])
#Ensure that this column is an integer so that it can be sorted correctly.
#Sort the dataframe based on fold change
foutdf = outdf.sort_values(by=["Support", "Fold Change Duplications"], ascending=False)

fname = OGsDir+"/FC-Dups_in_"+TaxClass+"_Orthogroups.txt"
#Write this dataframe to a \t seperated .txt    
foutdf.to_csv(fname, sep="\t", index=False)

ClassSpecOGfn = OGsDir+"/%s_Specific_OGs.txt" % (TaxClass)
ClassSpecLossfn = OGsDir+"/%s_Lost_OGs.txt" % (TaxClass)
with open(ClassSpecOGfn,'w') as ClassSpecOGf:
	for O in ClassSpecOGs:
		ClassSpecOGf.write("%s\n" % O)
with open(ClassSpecLossfn,'w') as ClassSpecLossf:
	for O in ClassSpecLoss:
		ClassSpecLossf.write("%s\n" % O)

#Do a heatmap of genes/orthogroup for each species from file /Orthogroups/Orthogroups.GeneCount.tsv
#Use subset of OGs from here.
#First get subset df. Then for heatmap, should be same as other script. 
#Don't use corr. Absolute numbers more interesting here.  

#Get average no. of genes/orthogroup for Class
OGdf = pd.read_csv(OG_count, sep="\t", index_col="Orthogroup")
OGdf.drop('Total', axis=1, inplace=True)
OGdf = OGdf[OGdf.index.isin(OGs)]
Classes = ['Bivalvia','Caudofoveata','Cephalapoda','Gastropoda','Monoplacophora','Polyplacophora','Scaphopoda','Solenogastres']
palette = sns.color_palette()
lut = dict(zip(map(str, Classes), palette))
# Convert the palette to vectors that will be drawn on the side of the matrix
# Mathes each class used to a colour. 
Spec_Classes = (OGdf.columns.to_series()).map(SpeClaD)
#Gets the index of the df (species names) and maps them to corresponding taxa Class
col_colors = Spec_Classes.map(lut)

ClusterGrid = sns.clustermap(OGdf, col_colors=col_colors)
#This creates the clustermap based on correlation
plt.setp(ClusterGrid.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
#This makes the y-axis titles horizontal
fig = plt.gcf()
filename = OGsDir+'/'+TaxClass+"_GeneDup_OGs_Heatmap.png"
fig.savefig(filename, dpi=300)
