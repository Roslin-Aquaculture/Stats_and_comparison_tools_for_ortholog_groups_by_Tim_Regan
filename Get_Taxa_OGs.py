#!/usr/bin/python
# Name/date: TimRegan/2022.08.11
# File: Get_Taxa_OGs.py
# Purpose: To ID and extract the names of orthogroups which have gene duplications retained in at least 50% of Taxa analysed. 
# Edit "Support threshold" for Duplications if 0.5 not desired.
# Also outputs taxa specific OGs, OGs whcih are in every species apart form this Taxa.
# Make sure to navigate to Orthofinder folder containing results folders.
# IMPORTANT: Edit as appropriate for your analyses and taxa on lines: 32, 41, 50, 61, 62, 69, 
# Usgae: python Get_Taxa_OGs_Lobster_v_Decapoda.py <Results_folder> <Node> <Taxa>
# e.g. 
# python Get_Taxa_OGs.py Results_RAXML-NG_Tree_1 N8 Lobster
# python Get_Taxa_OGs.py Results_RAXML-NG_Tree_1 N3 Decapoda
# python Get_Taxa_OGs.py Results_May16 N7 Bivalvia
# python Get_Taxa_OGs.py Results_May16 N8 Gastropoda

from __future__ import division
import matplotlib
matplotlib.use('Agg')
import scipy
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

# Dictionary to rename species names from file names, if necessary
File_species = { "D_mela_protein" :"Drosophila_melanogaster", \
"D_magn_protein" : "Daphnia_magna", \
"H_azte_protein" : "Hyalella_azteca", \
"H_amer_protein" : "Homarus_americanus", "H_gama_protein" : "Homarus_gammarus", \
"P_chin_protein" : "Penaeus_chinensis", "P_japo_protein" :"Penaeus_japonicus", "P_mono_protein" : "Penaeus_monodon", "P_vann_protein" : "Penaeus_vannamei", \
"P_trit_protein" : "Portunus_trituberculatus", \
"P_clar_protein" : "Procambarus_clarkii", "C_dest_protein" : "Cherax_destructor"}

# Taxa groups - state which species are in which taxonomic group following examples included here:
Fly = ["Drosophila_melanogaster"]
Water_flea = ["Daphnia_magna"]
Amphipod = ["Hyalella_azteca"]
Lobster = ["Homarus_americanus", "Homarus_gammarus"]
Shrimp = ["Penaeus_vannamei", "Penaeus_chinensis", "Penaeus_japonicus", "Penaeus_monodon"]
Crab = ["Portunus_trituberculatus"]
Crayfish = ["Procambarus_clarkii", "Cherax_destructor"]

# Create a dictionary using information provided in taxa groups following this example:
SpeClaD = { "Drosophila_melanogaster" : "Fly", \
"Daphnia_magna" : "Water_flea", \
"Hyalella_azteca" : "Amphipod", \
"Homarus_americanus" : "Lobster", "Homarus_gammarus" : "Lobster", \
"Penaeus_chinensis" : "Shrimp", "Penaeus_japonicus" : "Shrimp", "Penaeus_monodon" : "Shrimp", "Penaeus_vannamei" : "Shrimp", \
"Portunus_trituberculatus" : "Crab", \
"Procambarus_clarkii" : "Crayfish", "Cherax_destructor" : "Crayfish"}

# SpeClaD = {'Acanthochitona_crinita': 'Polyplacophora', 'Achatina_fulica': 'Gastropoda', 'Aplysia_californica': 'Gastropoda', 'Bathymodiolus_platifrons': 'Bivalvia', 'Biomphalaria_glabrata': 'Gastropoda', 'Crassostrea_gigas': 'Bivalvia', 'Crassostrea_virginica': 'Bivalvia', 'Cristaria_plicata': 'Bivalvia', 'Elysia_chlorotica': 'Gastropoda', 'Gadila_tolmiei': 'Scaphopoda', 'Gymnomenia_pellucida': 'Solenogastres', 'Laevipilina_hyalina': 'Monoplacophora', 'Laternula_elliptica': 'Bivalvia', 'Lottia_gigantea': 'Gastropoda', 'Lymnaea_stagnalis': 'Gastropoda', 'Mizuhopecten_yessoensis': 'Bivalvia', 'Modiolus_philippinarum': 'Bivalvia', 'Mya_arenaria': 'Bivalvia', 'Mya_truncata': 'Bivalvia', 'Mytilus_edulis': 'Bivalvia', 'Mytilus_galloprovincialis': 'Bivalvia', 'Octopoteuthis_deletron': 'Cephalapoda', 'Octopus_bimaculoides': 'Cephalapoda', 'Octopus_sinensis': 'Cephalapoda', 'Octopus_vulgaris': 'Cephalapoda', 'Pecten_maximus': 'Bivalvia', 'Pinctada_fucata': 'Bivalvia', 'Pomacea_canaliculata': 'Gastropoda', 'Scutopus_ventrolineatus': 'Caudofoveata', 'Vampyroteuthis_infernalis': 'Cephalapoda', 'Wirenia_argentea': 'Solenogastres'}

# State collective groups of taxa for analysing following the examples included here e.g.:
Decapoda = Lobster+Shrimp+Crab+Crayfish
All_taxa = Decapoda+Fly+Water_flea+Amphipod

#Mollusca = Caudofoveata+Bivalvia+Cephalopoda+Gastropoda+Monoplacophora+Polyplacophora+Scaphopoda+Solenogastres
#Classes = ['Bivalvia','Caudofoveata','Cephalapoda','Gastropoda','Monoplacophora','Polyplacophora','Scaphopoda','Solenogastres']

node = str(sys.argv[2])
# This is the  root node from the species tree. 
Support_threshold = 0.5

TaxGroup = str(sys.argv[3])
TaxList = eval(TaxGroup)
# NonTaxa = [s for s in Decapoda if s not in TaxList]
# For Decapod branch vs. other decapods i.e. conparing one tree branch to its surrounding branches without necessarily including all members
NonTaxa = [s for s in All_taxa if s not in TaxList]
# For any taxa group listed above versus everything else included in tree

dups = str(sys.argv[1])+'/Gene_Duplication_Events/Duplications.tsv'
#The Duplications.tsv file
#Load in Duplications.tsv filtered by node of interest
dupf = pd.read_csv(dups, sep="\t", index_col=None, header=0)
dupf = dupf.rename(columns=File_species)
#Renames columns according to species names
is_node =  dupf['Species Tree Node']==node
dupf_node = dupf[is_node]

OG_count = str(sys.argv[1])+'/Orthogroups/Orthogroups.GeneCount.tsv'
#Load in this file (counts per OG/species)
countdf = pd.read_csv(OG_count, sep="\t", index_col="Orthogroup", header=0)

countdf = countdf.rename(columns=File_species)
#Renames columns according to species names

#Get average no. of genes/orthogroup for Taxa group
Cldf = countdf.loc[:, countdf.columns.isin(TaxList)]
ClSpecDf = Cldf[(Cldf != 0).all(1)]
ClSpecDf = list(ClSpecDf.index.values)
ClSpecLoss = Cldf[(Cldf == 0).all(1)]
ClSpecLoss = list(ClSpecLoss.index.values)

NCdf = countdf.loc[:, countdf.columns.isin(NonTaxa)]
ClSpecNDf = ClSpecNDf = NCdf[(NCdf == 0).all(1)]
ClSpecNDf = list(ClSpecNDf.index.values)
NClSpecG = NCdf[(NCdf != 0).all(1)]
NClSpecG = list(NClSpecG.index.values)

TaxaSpecOGs = set(ClSpecDf).intersection(ClSpecNDf)
TaxaSpecOGs = list(sorted(TaxaSpecOGs))
TaxaSpecLoss = set(ClSpecLoss).intersection(NClSpecG)
TaxaSpecLoss = list(sorted(TaxaSpecLoss))

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
	avgTx = Cldf.loc[OG,'mean']
	avgnTx = NCdf.loc[OG,'mean']
	if avgnTx == 0:
		FC = float("inf")
	else:
		FC = avgTx/avgnTx
		FC = "{:.2f}".format(FC)
	return FC

def OG_MWU(OG):
    pvalue = None
    TxCounts = Cldf.loc[OG]
    nTxCounts = NCdf.loc[OG]
    implicit_count_1 = [count for count in TxCounts if count > 0]
    implicit_count_2 = [count for count in nTxCounts if count > 0]
    if len(implicit_count_1) >= 2 and len(implicit_count_2) >= 2:
        if len(set(implicit_count_1)) == 1 and len(set(implicit_count_2)) == 1 and set(implicit_count_1) == set(implicit_count_2): # equal
            pvalue = 1.0
        else:
            try:
                pvalue = scipy.stats.mannwhitneyu(implicit_count_1, implicit_count_2, alternative="two-sided")[1]
            except ValueError: # throws ValueError when all numbers are equal
                pvalue = 1.0
    return pvalue


#Get orthogroups with duplications common to more than half of Taxa ("OGs"). 
OGs = get_Node_OGs(Support_threshold)

#Make a corresponding list of corresponding fold change and max Support value of Bv gene duplications for each OG 
#For bonus, copy .fa for each orthogroup to new directory
OGsDir = 'OGs_'+str(sys.argv[1])
if not os.path.exists(OGsDir):
    os.makedirs(OGsDir)
TaxaOGsDupsDir = OGsDir + '/' + TaxGroup + "_Dups"
if not os.path.exists(TaxaOGsDupsDir):
    os.makedirs(TaxaOGsDupsDir)
FC = []
Supp = []
MWU =[]
for OG in OGs:
	FC.append(OG_FC(OG))
	Supp.append(OG_Supp(OG))
	MWU.append(OG_MWU(OG))
	fa = str(sys.argv[1])+'/Orthogroup_Sequences/'+OG+'.fa'
	fac = TaxaOGsDupsDir+'/'+OG+'.fa'
	copyfile(fa, fac)

TaxaOGsLossDir = OGsDir + '/' + TaxGroup + "_Loss"
if not os.path.exists(TaxaOGsLossDir):
    os.makedirs(TaxaOGsLossDir)
for O in TaxaSpecLoss:
	fa = str(sys.argv[1])+'/Orthogroup_Sequences/'+O+'.fa'
	fac = TaxaOGsLossDir+'/'+O+'.fa'
	copyfile(fa, fac)

TaxaOGsSpecDir = OGsDir + '/' + TaxGroup + "_Specific"
if not os.path.exists(TaxaOGsSpecDir):
    os.makedirs(TaxaOGsSpecDir)
for O in TaxaSpecOGs:
	fa = str(sys.argv[1])+'/Orthogroup_Sequences/'+O+'.fa'
	fac = TaxaOGsSpecDir+'/'+O+'.fa'
	copyfile(fa, fac)

#Make an output file with orthogroups and fold change of duplications in TaxGroup vs. Non-TaxGroup
outdf = pd.DataFrame(list(zip(OGs,FC,Supp,MWU)), columns=["Orthogroup","Fold Change Duplications","Support", "MWU pvalue"])
outdf['Fold Change Duplications'] = pd.to_numeric(outdf['Fold Change Duplications'])
#Ensure the 'Fold Change Duplications' column is an integer so that it can be sorted correctly.
#Sort the dataframe based on fold change
foutdf = outdf.sort_values(by=["MWU pvalue"], ascending=True)

fname = OGsDir+"/FC-Dups_in_"+TaxGroup+"_Orthogroups.txt"
#Write this dataframe to a \t seperated .txt    
foutdf.to_csv(fname, sep="\t", index=False)

TaxaSpecOGfn = OGsDir+"/%s_Specific_OGs.txt" % (TaxGroup)
TaxaSpecLossfn = OGsDir+"/%s_Lost_OGs.txt" % (TaxGroup)
with open(TaxaSpecOGfn,'w') as TaxaSpecOGf:
	for O in TaxaSpecOGs:
		TaxaSpecOGf.write("%s\n" % O)
with open(TaxaSpecLossfn,'w') as TaxaSpecLossf:
	for O in TaxaSpecLoss:
		TaxaSpecLossf.write("%s\n" % O)

