#!/usr/bin/python
#Name/date: TimRegan/2020.04.14
#File: GraphOrthoStats.py
#Purpose: To make graphs of the comparitive statistics between species following Orthofinder analysis.
#Use in dir above Results folder
#Usgae: python GraphOrthoStats.py <Results_folder> <Row to graph e.g. % genes assigned>
#e.g. python GraphOrthoStats.py Results_Apr11 "Percentage of genes in species-specific orthogroups"
#or
#python GraphOrthoStats.py Results_Apr11 "Percentage of genes in orthogroups"
#or
#python GraphOrthoStats.py Results_May01 "Percentage of unassigned genes"

import sys
import matplotlib
matplotlib.use('Agg')
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

Stats = str(sys.argv[1])+'/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv'
#The Comparitive stats per species .tsv file

Graphed = sys.argv[2]
#The row from file to be graphed

df = pd.read_csv(Stats, sep='\t', index_col=0, header=0)
df = df.loc[Graphed, :]
#Data extracted as a series.

plt.rcParams["figure.figsize"] = (6, 4)
#sets figure size
fig = sns.barplot(df.index, df.values, linewidth=1, edgecolor=".2")
#sets data to use, orders the data alphabetically (see above for species sorted) sets indvidual points to 0.5 size. Line width is the thickness of boxplot lines/borders.
fig = fig.set_title(Graphed)
plt.xticks(rotation=30, ha="right", fontsize=8)
#x axis labels rotated and made size 8. 
plt.xlabel("Species")
plt.ylabel("% Total Genes")
plt.tight_layout()

fig = plt.gcf()
filename = Graphed+'.png'
fig.savefig(filename, dpi=300)
