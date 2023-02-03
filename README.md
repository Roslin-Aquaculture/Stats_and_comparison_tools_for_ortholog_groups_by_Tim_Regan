# Stats and comparison tools for ortholog groups !!!!!
Analysis of diversity across orthologous gene groups (orthogroups) between defined taxa (e.g. Bivalves vs. other Molluscs or Ostrea edulis vs. other Bivalves). 
To run this analysis, we first take a directory of peptide sequences.

Table of Contents
=================
* [Preparing Peptide sequences](#preparing-peptide-sequences)
* [Annotation of peptide assemblies](#annotation-of-peptide-assemblies)
* [Orthofinder analyses](#orthofinder-analyses)
* [Analysing OrthoFinder output](#analysing-orthofinder-output)
* [Taxa specific Orthogroups (OGs)](#taxa-specific-orthogroups-ogs)
* [Kinfin analyses](#kinfin-analyses)
* [Searching for proteins of interest](#searching-for-proteins-of-interest)
* [Synteny mapping of selected OG genes](#synteny-mapping-of-selected-og-genes)


# Preparing Peptide sequences
See: [Preparing_peptides.sh](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/Preparing_peptides.sh) \
Peptides are filtered for the longest isoform using [AGAT, available on Github](https://github.com/NBISweden/AGAT/blob/master/bin/agat_sp_keep_longest_isoform.pl).
Name files in the format "Genus_species.fa" (OrthoFinder uses filename without the extension as species name). \
Also rename fasta headers to include genus and species information - this avoids duplicate sequence name errors and allows easier species ID when examining genes later:\
`$awk '{ gsub(">",">G.species_"); print $0 }' Protein_sequences.fa > Genus_species.fa`

CD-Hit is then used to examine clusters of sequences (may also use this to filter for longest cluster sequence if using AGAT (agat_sp_keep_longest_isoform.pl) from the previous step is not possible).\
See: [cdhit_cluster_analysis.py](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/cdhit_cluster_analysis.py) \
These results can be visualised using [cdhit_cluster_plot.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/cdhit_cluster_plot.py)

BUSCO can also be used to assess and report completeness of assemblies. These are useful metrics when further analysing your orthogroup results e.g. prior knowledge of a large proportion of dulpicatied genes from BUSCO may explain a lot of supposed duplication events in ortholog analyses. 

# Annotation of peptide assemblies
Use interproscan to annotate all peptides used. 
A useful approach is to first concatenate all peptide sequences together, then split this file into many smaller files and run individually on Interproscan before finally combining the resulting output.\
**NOTE:** *Make sure that you are using the filtered and renamed peptide filesto avoid naming mismatch during later analysis* \
If you are working on teh University of Edinburgh computing cluster: See [Interproscan_array.sh](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/Interproscan_array.sh)\
Combine results into one file e.g. \
`cat *.tsv > All_proteins.tsv` \
This will be used later in Kinfin analysis

# Orthofinder analyses
See these tutorials for using OrthoFinder: https://davidemms.github.io/menu/tutorials.html

Run OrthoFinder according to tutorial and requirements. 

**Note:** *MSA is more appropriate to run, but may take much longer. If you have many assemblies (e.g. >8) Running other tree alignmers such as RAXML or IQTREE may be more appropriate for species tree, but not required for gene trees and will take an impossible length of time for most studies. Suggest generating species tree seperately if required.* 

# Analysing OrthoFinder output
To get an initial visualisation at what the overall results look like, use: [GraphOrthoStats.py](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/GraphOrthoStats.py)

# Taxa specific Orthogroups (OGs)
In this case, we wished to compare orthogroups/gene families across taxa e.g. Bivalvia vs. Mollusca or Ostrea edulis vs Bivalvia or Lobsters vs Decapoda. 

See: [Get_Taxa_OGs.py](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/Get_Taxa_OGs.py) \
Following instructions in script notes very carefully. 
Need to give some information to first Python file e.g. species within taxa group etc. 

To get Orthogroups (OGs) associated with a particular taxa group, run: \
`$python Get_Taxa_OGs.py <Results_folder> <Node> <Taxa>`

Where `<Results_folder>` is the Orthofinder results folder 
(e.g. Results_May22), \
`<Taxa>` is the Taxa group of interest (e.g. Bivalvia) \
and `<Node>` is the Species_tree_duplications node associated with your taxa.

**NOTE:** Make sure to define your specific taxa groupings with lists of individual taxa/species peptide sequences with in-line editing of the script *(see line 8 of script)*.  

This outputs all of the OGs specific to that taxa, lost by that taxa group, and OGs with gene duplications retained by >50% (default, can chnage under "support" value in script, anywhere from 0.0 - 1.0) of species within that Taxa.

From the OGs with duplication events, it produces a table with OGs ranked in order of how many taxa group members have the duplication (Support), and the avg. Fold-Change of genes within the OG compared to other taxa in addition to a Mann-Whitney U test to provide a p value for OGs significantly expanded in your taxa of interest vs other taxa included in the comparison e.g. labster vs. other Decapods or Bivalves vs. other Molluscs. 

# Kinfin analyses 
KinFin is a toolkit for the local implementation of automated analyses of clustered proteins through custom taxon grouping:

it takes a protein clustering output by tools such as OrthoFinder or OrthoMCL, alongside functional annotation data (InterProScan output), and user-defined species taxonomies, to derive rich aggregative annotation of orthology groups.

Please see https://kinfin.readme.io/docs for more information. 

**Note:** takes in interproscan annotation file: \
 `iprs2table.py -i Clustered_proteins_IPR.tsv --domain_sources GO,Pfam`

 An example of code to run Kinfin (requires manually created config file, see: https://kinfin.readme.io/docs/starting-from-a-clustering-files) \
 
 ./kinfin/kinfin \  
  -g OrthoFinder/Results_Aug04/Orthogroups/Orthogroups.txt \  
-c config.txt \
-s OrthoFinder/Results_Aug04/WorkingDirectory/SequenceIDs.txt \
-p OrthoFinder/Results_Aug04/WorkingDirectory/SpeciesIDs.txt\
-f functional_annotation.txt \
-a /protein_files \
-t Rooted_raxml-ng_tree.txt \
--infer_singletons -o Kinfin_default

# Searching for proteins of interest
We might be looking for which orthgroups certain proteins have fallen into (or which other proteins they are in an orthogroup with). e.g. TLR proteins, IFN genes etc. \
See [Find_containing.py](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/Find_containing.py) \
Requires a directory of .fa files to search for protein description.

# Synteny mapping of selected OG genes
Synteny mapping can be performed to compare Orthogroups between species. 
Best done using species which have good genome assemblies ideally at a chromosomal level. 
To obtain Chromosomal location of genes of interest for Synteny mapping, see: [Get_OG_member_loc.py](https://github.com/Roslin-Aquaculture/Stats_and_comparison_tools_for_ortholog_groups_by_Tim_Regan/blob/master/Get_OG_member_loc.py)

*Requires the dictionary of unique identifiers used previously for each species header (see above). 
It outputs two results files for each species of interest: the gene names and their sequences specific to each species.
The list of gene names can then be run against a .gff annotaion file to give the chromosome, gene start, gene end and gene name. See notes at end of script. 
