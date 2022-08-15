# Orthofinder_Pipeline
Orthogroup and gene group diversification/selection analysis between groups of species (e.g. Bivalves vs. other Molluscs or Ostrea edulis vs. other Bivalves). 
To run this analysis, we first take a directory of peptide sequences.

# Preparing Peptide sequences
See: [Preparing_peptides.sh](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/Preparing_peptides.sh) \
Peptides are filtered for the longest isoform using [AGAT, available on Github](https://github.com/NBISweden/AGAT/blob/master/bin/agat_sp_keep_longest_isoform.pl).
Name files in the format "Genus_species.fa" (OrthoFinder uses filename without the extension as species name). \
Also rename fasta headers to include genus and species information - avoids duplicate sequence name errors:\
`$awk '{ gsub(">",">G.species_"); print $0 }' Genus_species.fa > Genus_species_renamed.fa`

CD-Hit is then run to examine clusters of sequences (may also use this to filter for longest cluster sequence if using AGAT is not possible).\
See: [cdhit_cluster_analysis.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/cdhit_cluster_analysis.py) \
These results can be plotted using [cdhit_cluster_plot.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/cdhit_cluster_plot.py)

BUSCO may also be run to assess completeness of assemblies. 

# Annotation of peptide assemblies
Use interproscan to annotate all peptides used. 
A useful approach is to first concatenate all (filtered, header renamed) peptide sequences together, then split this file into many smaller files and run individually on Interproscan before finally combining the resulting output.\
See [Interproscan_array.sh](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/Interproscan_array.sh)\
Combine results into one file e.g. \
`cat *.tsv > All_proteins.tsv` \
This will be used later in Kinfin analysis

# Orthofinder analyses
See these tutorials for using OrthoFinder: https://davidemms.github.io/menu/tutorials.html

Run OrthoFinder according to tutorial and requirements. 

Note: MSA is more appropriate to run, but may take longer. Running other tree alignmers such as RAXML or IQTREE may be more appropriate for species tree, but not required for gene trees and will take an impossible length of time for most studies. Suggest generating species tree seperately if required. 

# Analysing OrthoFinder output
To get an initial look at what the results look like, see: [GraphOrthoStats.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/GraphOrthoStats.py)

# Class specific Orthogroups (OGs)
In this case, we wished to compare orthogroups/gene families across taxa e.g. Bivalvia vs. Mollusca or Ostrea edulis vs Bivalvia or Lobsters vs Decapoda. 

See: [Get_Class_OGs.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/Get_Class_OGs_Lobster_v_Decapoda.py) \
Following instructions in script notes very carefully. 
Need to give some information to first Python file e.g. species within class group etc. 

To get Orthogroups (OGs) associated with a particular Class, run: \
`$python Get_Class_OGs.py <Results_folder> <Node> <Class>`

Where `<Results_folder>` is the Orthofinder results folder 
(e.g. Results_May22), \
`<Class>` is the Class of interest (e.g. Bivalvia) \
and `<Node>` is the Species_tree_duplications node associated with your Class.

This outputs all of the OGs specific to that class, lost by that class, and OGs with gene duplications retained by >50% (default, can chnage under "support" value in script, anywhere from 0.0 - 1.0) of species within that Class.

From the OGs with duplication events, it produces a table with OGs ranked in order of how many Class members have the duplication (Support), and the avg. Fold-Change of genes within the OG compared to other Classes in addition to a Mann-Whitney U test to provide a p value for OGs significantly expanded in your taxa of interest vs other taxa included in the comparison e.g. labster vs. other Decapods or Bivalves vs. other Molluscs. 

# Kinfin analyses 
KinFin is a toolkit for the local implementation of automated analyses of clustered proteins through custom taxon grouping:

it takes a protein clustering output by tools such as OrthoFinder or OrthoMCL, alongside functional annotation data (InterProScan output), and user-defined species taxonomies, to derive rich aggregative annotation of orthology groups.

Please see https://kinfin.readme.io/docs for more information. 

Note: takes in interproscan annotation file: \
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
See [Find_containing.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/Find_containing.py) \
Requires a directory of .fa files to search for protein description.

# Making a Network Graph of OG genes
The genes trees provided by OrthoFinder are very useful for identifying Orthologs vs. Paralogs etc.
However, these can become confusing to the point of being unreadable when multiple species and genes are present. \
To view the OG using a Network Graph, see: [OGfa2graph.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/OGfa2graph.py)

Runs a BLAST all v. all on each of the OG sequences and creates a pairwise.txt file of their relationships. \
Requires taxonomic species/grouping dictionary from before for colouring. \
Also requires a dictionary stating the unique header identifiers from each species assembly
i.e. \
Headers = {'ID1':'Species1', 'IDn':'Speciesn'} \
e.g. Headers = {'Cgiga':'Crassostrea_gigas', 'Hgamm':'Homarus_gammarus'} etc. 

Also requires Blast to be callable from the home directory. Make sure the module is loaded, or else edit the script. 

Options for stringency of comparisons are available. 
This way, we can delve into when gene duplications occurred (phyogenetic relationships), which are Species/Class specific.
Colour coding the graph by species/class makes this much more comprehensible and the graph can be filtered according to similarity.nearest neighbour etc.

# Synteny mapping of selected OG genes
Synteny mapping can be performed to compare Orthogroups between species. 
Best done using species which have good genome assemblies ideally at a chromosomal level. 
To obtain Chromosomal location of genes of interest for Synteny mapping, see: [Get_OG_member_loc.py](https://github.com/Roslin-Aquaculture/Orthofinder_Bivalvia/blob/master/Get_OG_member_loc.py)

Requires the dictionary of unique identifiers used previously for each species header (see above). 
It outputs two results files for each species of interest: the gene names and their sequences specific to each species.
The list of gene names can then be run against a .gff annotaion file to give the chromosome, gene start, gene end and gene name. See notes at end of code. 




