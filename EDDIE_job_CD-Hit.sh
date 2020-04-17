#!/bin/bash
#$ -N TR_cdhit_1
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=10G
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 24 hrs: -l h_rt
#  memory limit of 10 Gbyte: -l h_vmem
#Tim_Regan
#2020.04.02
#Compares number of peptide sequences in assembly to number of clusters calculated using CD-HIT
#For more, see https://bioinformaticsreview.com/20190922/how-to-cluster-peptide-protein-sequences-using-cd-hit-software/

# Initialise the environment modules
. /etc/profile.d/modules.sh

#load required modules
module load igmm/apps/cdhit/4.6.8

#go to peps directory
cd /exports/eddie/scratch/tregan/Mollusc_peps
mkdir cd-hit

# run cd-hit
for f in *.fa ; do cd-hit -i "$f" -o cd-hit/"$f".cdhit -T 4 ; done

#Make a file with species name on each line
ls *.fa | awk '{print substr($0,1,length()-3);}' > Species.txt

#for line in "Species.txt", get number of peptide seqs and clusters
while read line
do
	grep -c "^>" "$line.fa" >> peps.txt
	grep -c "^>" "cd-hit/$line"*.clstr >> clstrs.txt
done < Species.txt

#for 1st column in "Species.txt" ($1), get number of clusters in $1*.clstr and add in column 3
paste Species.txt peps.txt clstrs.txt | column -s $'\t' -t > Species.peps.clstrs.txt

#Finally produce a file containing columns for each species, 1. Species, 2. No. of peptide seqs, 3. No. of CD-Hit clusters, and 4. No. of seqs/Cluster.
awk '{ print $1, "\t", $2 / $3 }' Species.peps.clstrs.txt > cd-hit.txt


rm Species.peps.clstrs.txt clstrs.txt peps.txt Species.txt

#For more analyses and graphs, see python scripts "cdhit_cluster_analysis.py" and "cdhit_cluster_plot.py"
