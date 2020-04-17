#!/bin/sh
#$ -N TR_kofamscan_2
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=10G
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 24 hrs: -l h_rt
#  memory limit of 10 Gbyte: -l h_vmem
# Initialise the environment modules
. /etc/profile.d/modules.sh

cd /exports/eddie3_homes_local/tregan

# load anaconda and activate
module load anaconda
source activate /exports/eddie/scratch/tregan/Kofamscan/kofamscan

cd /exports/eddie/scratch/tregan/Kofamscan

# run kofamscan
#usage: ./exec_annotation -o result.txt query.fasta
#this outputs the results file to the current location. 

cd /exports/eddie/scratch/tregan/Kofamscan
for f in /exports/eddie/scratch/tregan/Mollusc_peps/OrthoFinder/BivalveOGs/*.fa ; 
do kofamscan/bin/exec_annotation -o "$f".txt $f ; done

tail -n +3 /exports/eddie/scratch/tregan/Mollusc_peps/OrthoFinder/BivalveOGs/*.fa.txt | cat > /exports/eddie/scratch/tregan/Mollusc_peps/OrthoFinder/BivalveOGs/N9_OGs_kofamscan.txt
#this skips the top 2 lines of every file and returns a concatenated file of all outputs (no headers).

# qsub -pe sharedmem 4 -R y /exports/eddie3_homes_local/tregan/JobScripts/TR_kofamscan_2.sh
