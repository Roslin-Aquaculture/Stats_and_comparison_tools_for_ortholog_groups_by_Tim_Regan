#!/bin/sh
#$ -N TR_busco_1
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

#NB: Ensure paths are stated in the config.ini file within BUSCO, ensure datasets are downloaded e.g. Mollusca/Metazoa etc.
#Set paths in config.ini:
#[augustus]
# path to augustus
#path = /exports/eddie3_apps_local/apps/community/roslin/augustus/3.3.32/bin/
#[etraining]
# path to augustus etraining
#path = /exports/eddie3_apps_local/apps/community/roslin/augustus/3.3.32/bin/
# path to augustus perl scripts, redeclare it for each new script
#[gff2gbSmallDNA.pl]
#path = /exports/eddie3_apps_local/apps/community/roslin/augustus/3.3.32/scripts
#[new_species.pl]
#path = /exports/eddie3_apps_local/apps/community/roslin/augustus/3.3.32/scripts
#[optimize_augustus.pl]
#path = /exports/eddie3_apps_local/apps/community/roslin/augustus/3.3.32/scripts
#[hmmsearch]
# path to HMMsearch executable
#path = /exports/eddie3_apps_local/apps/community/roslin/hmmer/3.1b2/binaries/

# load anaconda and activate
module load anaconda
source activate busco

#load required modules
module load roslin/hmmer/3.1b2 
module load roslin/blast+/2.9.0
module load roslin/augustus/3.3.32
module load roslin/R/3.4.3

#go to busco directory
cd /exports/eddie/scratch/tregan/busco/busco
cd /exports/eddie/scratch/tregan/Mollusc_peps/BUSCO

# run busco
for f in *.fa ; do python ../../busco/busco/scripts/run_BUSCO.py -m proteins -i "$f"  -o "$f".out -l ../../busco/mollusca_odb10 ; done
python ../busco/busco/scripts/run_BUSCO.py -m proteins -i Crassostrea_gigas.v9.fa  -o Crassostrea_gigas.v9.ml.out -l ../busco/mollusca_odb10

