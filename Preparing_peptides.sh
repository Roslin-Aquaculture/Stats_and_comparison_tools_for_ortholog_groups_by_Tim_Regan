#!/bin/bash
#Tim_Regan
#2022.04.01
#Preparing_peptides.sh
#Clusters peptide sequences for each species and retains the longest member of each cluster as a proxy for longest isoform of each gene. 

#Prior to running this, first ensured headers all contained first letter of genus and first 4 letters of species, e.g.
#$ awk '{ gsub(">",">Medul_"); print $0 }' Mytilus_edulis.fa  > Mytilus_edulis_.fa
#Repeat for each species

#load required modules
module load igmm/apps/cdhit/4.6.8

#go to peps directory
cd /exports/eddie/scratch/tregan/Mytilus
mkdir cd-hit


# run cd-hit on each peptome
for peps in *.fa; do  cd-hit -i $peps -o ../"$peps" -T 4; done

mkdir Ortho 
mv *.fa Ortho/

cd Ortho

#Make a file with species name on each line
ls *.fa | awk '{print substr($0,1,length()-3);}' > Species.txt

#for line in "Species.txt", get number of peptide seqs and clusters
while read line; do; grep -c "^>" "$line.fa" >> peps.txt; grep -c "^>" "cd-hit/$line"*.clstr >> clstrs.txt; done < Species.txt

#for 1st column in "Species.txt" ($1), get number of clusters in $1*.clstr and add in column 3
paste Species.txt peps.txt clstrs.txt | column -s $'\t' -t > Species.peps.clstrs.txt

awk '{ print $1, "\t", $2 / $3 }' Species.peps.clstrs.txt > cd-hit.txt
rm Species.peps.clstrs.txt clstrs.txt peps.txt Species.txt


#Genes per assembly
grep -c '>' *.fa >> Counts.txt
#e.g.
#$ cat Counts.txt
#Bathymodiolus_platifrons.fa:29260
#Crassostrea_gigas.fa:27156
#Crassostrea_virginica.fa:26306
#Cristaria_plicata.fa:29702
#Laternula_elliptica.fa:34446
#Mizuhopecten_yessoensis.fa:22816
#Modiolus_philippinarum.fa:32873
#Mya_arenaria.fa:24318
#Mya_truncata.fa:44385
#Mytilus_coruscus.fa:52829
#Mytilus_edulis.fa:34746
#Mytilus_galloprovincialis.fa:37783
#Ostrea_edulis.fa:31739
#Pecten_maximus.fa:23651
#Pinctada_fucata.fa:29892
#Scapharca_broughtonii.fa:22451

