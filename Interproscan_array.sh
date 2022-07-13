#!/bin/sh
#$ -N IPR_array
#$ -cwd
#$ -l h_rt=48:00:00
#$ -l h_vmem=64G
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 48 hrs: -l h_rt
#  memory limit of 64 Gbyte: -l h_vmem
# Initialise the environment modules
. /etc/profile.d/modules.sh
module load roslin/interproscan/5.33-72.0 

cd /exports/eddie/scratch/tregan/IPR_data

F=`sed -n ${SGE_TASK_ID}p < /exports/eddie/scratch/tregan/IPR_data/IPR_List.txt`
#Where "IPR_List.txt" contains a list of each protein file on a new line. 

#Create a new output folder and temp folder for each file:
mkdir `basename $F .fa`
mkdir Temp_"$SGE_TASK_ID"

#Exectute interproscan command:
interproscan -i $F \
-d "$F"/ \
--tempdir Temp_"$SGE_TASK_ID"/ \
-cpu 16 \
-dp \
-t p \
--goterms \
-appl TIGRFAM,ProSiteProfiles,Pfam \
-f TSV

rm -r Temp_"$SGE_TASK_ID"

#Submit job with number of jobs to be run at once, e.g. if there are 34 individual protein files, use the following:
#qsub -t 1-34 IPR_array.sh
