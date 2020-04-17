#!/bin/sh
#$ -N TR_orthofinder_1
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=10G
#  These options are:
#  job name: -N
#  use the current working directory: -cwd
#  runtime limit of 24 hrs: -l h_rt
#  memory limit of 10 Gbyte: -l h_vmem

#no-arch v2.3.11 available here: https://anaconda.org/bioconda/orthofinder/files
#May need to manually install no-arch version.

# Initialise the environment modules
. /etc/profile.d/modules.sh
cd /exports/eddie3_homes_local/tregan
# load anaconda and activate

module load anaconda
source activate Ortho

#Go to peptides dir and run orthofinder
cd /exports/eddie/scratch/tregan
orthofinder -S diamond -f Mollusc_peps/ -t 8
