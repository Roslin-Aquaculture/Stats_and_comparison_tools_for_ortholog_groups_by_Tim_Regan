#!/usr/bin/python
#Name/date: TimRegan/2020.04.28
#File: Get_OG_member_loc.py
#Purpose: This script is used to find the location within a genome assembly of each genes within an orthogroup  
#Usgae: python Get_OG_member_loc.py <infile>
#Note: ensure that line 27 is updates to reflect fasta headers relating to each species you are using. Some examples are provided.

"""
1. Select OGs and corresponding Gene Tree. 
2. Retrieve species specific sequences/gene names from selected OGs
3. Map these genes to the genomes(tblastn)/find chr location from .gff:
4. Use these coordinates along with Chr length for each species to create CIRCOS plots

Again, each seq is labeled with Species ID.E.g.: 
awk '{ gsub(">",">G.species_"); print $0 }' Genus_species.fa > Genus_species_renamed.fa
<for $file in Mollusc_peps/*.fa ; do head -n 1 $file ; done> 
to see uniq seq headers. Add these to dict here. Use in later analyses. 

Some assemblies have starts_with "G_species_"
Other assemblies (e.g. MDB) have unique identifiers at start of gene name. 
"""
import sys
from Bio import SeqIO

infile = sys.argv[1]

Headers = {"Cgiga_":"C.gigas", "Cvirg_":"C.virginica", "Pmaxi_":"P.maximus", "Pcana_":"P.canaliculata", "Ovulg_":"O.vulgaris"}
#A dictionary containing the associated patterns specific to each species in the fasta headers. 

def return_seqs(species, pattern, fa):
	"""Requires dictionary with species specific header patterns in each sequence of fasta file
	returns a fasta file containing all sequences specific to a given species from a given fasta"""
	outseq = "%s_%s" % (species, fa)
	outnames = outseq.replace(".fa", "_genes.txt")
	with open(outnames, 'a') as outgenes:
		with open(outseq, 'a') as outfa:
			for seq_record in SeqIO.parse(fa, "fasta"):
				if pattern in seq_record.description:
					specseq = seq_record.format("fasta")
					specgen = "%s\n" % (seq_record.id)
					specgen =specgen.strip(pattern)
					outfa.write(specseq)
					outgenes.write(specgen)

for pattern in Headers:
	species = Headers[pattern]
	return_seqs(species, pattern, infile)

"""
This works nicely. Delivers gene names and seqs. 
The gene names can then be run against a GFF to ID location. 
#for genenames in *genes.txt; do grep -f $genenames < *.gff | cut -f1,4,5,9 | cut -d";" -f1 | awk -F'ID=cds-' '{print $1 $2}' > Loc."$genenames"; done
#takes a file, gene_names.txt, passes each line (gene) through the gtf annotation file and gives Chr loc.

Retrieve .fna, .faa and .gff for each species:
e.g.
Pecten maximus	xPecMax1.1	GCA_902652985.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/6579/100/GCF_902652985.1_xPecMax1.1/GCF_902652985.1_xPecMax1.1_genomic.gff.gz
https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/6579/100/GCF_902652985.1_xPecMax1.1/GCF_902652985.1_xPecMax1.1_genomic.gtf.gz
Crassostrea virginica	C_virginica-3.0	GCF_002022765.2	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/022/765/GCF_002022765.2_C_virginica-3.0/GCF_002022765.2_C_virginica-3.0_genomic.gff.gz
Pomacea canaliculata	ASM479433v1	GCA_004794335.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/073/045/GCF_003073045.1_ASM307304v1/GCF_003073045.1_ASM307304v1_genomic.gff.gz
Octopus sinensis/vulgaris	ASM634580v1	GCF_006345805.1	https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/345/805/GCF_006345805.1_ASM634580v1/GCF_006345805.1_ASM634580v1_genomic.gff.gz
"""

"""This works nicely. Delivers gene names and seqs. 
The gene names can then be run against a GFF to ID location:
#for genenames in *genes.txt; do grep -f $genenames < *.gff | cut -f1,4,5,9 | cut -d";" -f1 | awk -F'ID=cds-' '{print $1 $2}' > Loc."$genenames"; done
takes a file, gene_names.txt, passes each line (gene) through the gtf annotation file and gives Chr loc.
Use std.out or something to put this into command line: sys.stdout.write

Use locations (start, stop, end, chr etc.) for Circos: 
https://training.galaxyproject.org/training-material/topics/visualisation/tutorials/circos/tutorial.html#introduction
https://usegalaxy.org/

"""
