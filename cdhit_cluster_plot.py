#!/usr/bin/python
#Name/date: TimRegan/2020.04.07
#File: cdhit_cluster_plot.py
#Purpose: To output the number of clusters, and sequences per cluster for a .cdhit.clstr cd-hit output file 
#Usage: python cdhit_cluster_plot.py

import sys
import matplotlib
matplotlib.use('Agg')
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

cdhitfs = input("Enter your directory containing seqs/cluster: ")
path = r''+cdhitfs # use your path
all_files = glob.glob(path + "/*.csv")

frame = pd.DataFrame()

for filename in all_files:
	df = pd.read_csv(filename, sep='\t', index_col=None, header=0)
	df = df.iloc[:, 1:2]
	#This just selects the column with species name and seqs/scluster (not cluster names)
	frame = pd.concat([frame, df], axis=1)
	#axis=1 tells concat to add columns side by side.
frame.to_csv("All_species_clstrs.csv", sep='\t')

#frame is a dataframe with 22 columns, headers named after species, and 77493 rows, all clusters.
#Get list of species from column headers and sort alphabetically

frame = pd.read_csv("All_species_clstrs.csv", sep='\t', index_col=0, header=0)
species = frame.columns.values.tolist()
speciesrt = sorted(species)

frame = pd.melt(frame)
#melts dataframe into 2 columns. This allows us to load teh entire dataset into seaborn - values vs. species. 

plt.rcParams["figure.figsize"] = (10, 4)
#sets figure size
fig = sns.stripplot(x="variable", y="value", data=frame, order=speciesrt, size=1, linewidth=0.1, jitter=0.4)
##sets data to use, orders the data alphabetically (see above for species sorted) sets indvidual points to size 1. Line width is the thickness of grey lines/borders.
fig = fig.set_title('Sequences/Cluster from CD-HIT analysis of Peptides')
plt.xticks(rotation=30, ha="right", fontsize=8)
#x axis labels rotated and made size 8. 
plt.xlabel("Species")
plt.ylabel("Sequences/Cluster")
plt.tight_layout()
fig = plt.gcf()
fig.savefig('cdhit_seqs.per.clstr.strip.png', dpi=300)

#The above code creates a strip plot. See code below for other viz options e.g. violin plots, boxplots etc. 
#Uncomment code accordingly. 

#plt.ylim(0, 6)
#This sets y axis limits to allow easier viz.
#fig1a = plt.gcf()
#fig1a.savefig('cdhit_seqs.per.clstr.box.cut.png', dpi=300)

#plt.rcParams["figure.figsize"] = (6, 4)
##sets figure size
#fig1 = sns.boxplot(x="variable", y="value", data=frame, order=speciesrt, fliersize=0.5, linewidth=0.1)
##sets data to use, orders the data alphabetically (see above for species sorted) sets indvidual points to 0.5 size. Line width is the thickness of boxplot lines/borders.
#fig1 = fig1.set_title('Sequences/Cluster from CD-HIT analysis of Peptides box')
#plt.xticks(rotation=30, ha="right", fontsize=8)
##x axis labels rotated and made size 8. 
#plt.xlabel("Species")
#plt.ylabel("Sequences/Cluster")
#plt.tight_layout()

#fig1 = plt.gcf()
#fig1.savefig('cdhit_seqs.per.clstr.box.png', dpi=300)

#plt.ylim(0, 6)
#This sets y axis limits to allow easier viz.
#fig1a = plt.gcf()
#fig1a.savefig('cdhit_seqs.per.clstr.box.cut.png', dpi=300)

### WARNING!!! ###
#Violin plots cannot be made at the same time as box plots here - pick one. Will fix code to do both at once eventually...

#plt.rcParams["figure.figsize"] = (6, 4)
#fig2 = sns.violinplot(x="variable", y="value", data=frame, order=speciesrt, cut=1, scale="count", width=1.5, linewidth=0.1)
#fig2 = fig2.set_title('Sequences/Cluster from CD-HIT analysis of Peptides')
#plt.xticks(rotation=30, ha="right", fontsize=8)
#plt.xlabel("Species")
#plt.ylabel("Sequences/Cluster")
#plt.tight_layout()

#fig2 = plt.gcf()
#fig2.savefig('cdhit_seqs.per.clstr.vln.png', dpi=300)

#plt.ylim(0, 8)
#fig2a = plt.gcf()
#fig2a.savefig('cdhit_seqs.per.clstr.vln.cut.png', dpi=300)


