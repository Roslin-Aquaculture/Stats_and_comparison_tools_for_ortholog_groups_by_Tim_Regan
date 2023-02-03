#!/usr/bin/python
#Name/date: TimRegan/2020.04.07
#File: cdhit_cluster_analysis.py 
#Purpose: To output the number of clusters, and sequences per cluster for a .cdhit.clstr cd-hit output file
#Usgae:	python cdhit_cluster_analysis.py <file>
#Or:	for f in *.cdhit.clstr ; do python cdhit_cluster_analysis.py $f ; done

import sys

def get_species(g):
	""" returns the species name from a file in format some_species.fa.cdhit.clstr """
	g = (str(g)).replace(".fa.cdhit.clstr","")
	return g


def count_clstrs(f):
	""" returns the number of clusters in a cd-hit .clstr output file """
	f = str(f)
	with open(f) as file:
		n = 0
		for line in file:
		    if line.startswith(">Cluster"):
		        n += 1
	return n

def count_seqs(f, s, c):
	""" takes a .cdhit.clstr file and outputs the sequences per cluster to a new file """
	""" where f = .cdhit.clstr file, g = species and h = no. of clusters in f  """
	counts = {}
	current = None
	with open(f) as fo:
	   for line in fo:
	       if line.startswith('>'):
	           current = line.strip()
	           counts[current] = 0
	       else:
	           counts[current] += 1
	#this assigns a count of sequences to each cluster.
	oot = s + ".csv"
	with open(oot,'a') as ot:
		topline = "%s\t%s\n" % (c, s)
		ot.write(topline)
		for entry, count in counts.items():
			entry = entry.replace(">Cluster ", "")
			ot.write('{} \t {:2d} \n'.format(entry, count))
			#This outputs the Cluster with the count to 2 decimal places. 

file = sys.argv[1]
species = get_species(file)
clusters = count_clstrs(file)
count_seqs(file, species, clusters)


#Next, graph seqs/cluster (box whiskers, violin plot and swarm/strip plots)
#Also graph clusters/seqs with no. of seqs as heading (column w/ over labeling)
#See "cdhit_cluster_plot.py"
