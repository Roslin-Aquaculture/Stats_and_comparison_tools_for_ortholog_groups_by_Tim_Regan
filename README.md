# Orthofinder_Bivalvia
Orthogroup and gene group diversification/selection analysis between groups of species, in this case, Bivalvia within Mollusca. 
Peptide assemblies analysed using CD-Hit and BUSCO prior to OrthoFInder analysis. 
Also, see these tutorials for using OrthoFinder: https://davidemms.github.io/menu/tutorials.html

To run this analysis, we first take a directory of peptide sequences.
Name files in the format "Genus_species.fa" (OrthoFinder uses filename without the extension as species name).
Run OrthoFinder. 

To get an initial look at what teh results look like, run
$python GraphOrthoStats.py <Results_folder> <Row to graph> 
Where <Results_folder> is the Orthofinder results folder (e.g. Results_May22)
<Row to graph> is e.g. "Percentage of genes in species-specific orthogroups" or "Percentage of genes in orthogroups"

In this case, we wished to compare Bivalvia within Mollusca i.e. a Class within a Phylum. 
Need to give some information to first Python file. 
Initially need to set each class used i.e. "Classes = ['Class1', 'Class2', 'Classn']"
Set each Class member i.e. "Class1 = ['Species1', 'Species2', 'Classn']"
Repeat for each Class. Then add these lists of classes together: "Phylum = Class1+Class2+Classn"

Finally, we need a dictionary telling the script which Class each Species is. 
"SpecClaD = {'Species1':'Class1', 'Species2':'Class1', 'Species3':'Class2', 'Speciesn':'Classn'}"

To get Orthogroups (OGs) associated with a particular Class, run:
$python Get_Class_OGs.py <Results_folder> <Node> <Class>
Where <Results_folder> is the Orthofinder results folder (e.g. Results_May22), 
<Class> is the Class of interest (e.g. Bivalvia)
<Node> is the Species_tree_duplications node associated with your Class.
This outputs all of the OGs specific to that class, lost by that class, and OGs with gene duplications retained by >50% of species within that Class.
It also gives a heatmap of the number of genes/species within each of the OGs with gene duplications in 50% of species within the Class. 
From these OGs with duplications, it produces a table with OGs ranked in order of how many Class members have the duplication (Support), and the avg. Fold-Change of genes within the OG compared to other Classes. 

To assign function to each of these Orthogroups, a peptide KEGG mapper called Kofamscan is used. 
This is the command line version of BlastKOALA. 
Format for file output:
Summarise these findings for each OG.
Make table of OGs from the OG_Duplications table (Ranked).

We can infer from this table the OGs which are most selected by or important to this Class. 
Looking at Class specific OGs or OGs lost by that Class are not ranked. 

At this point, you will have a good idea on which OGs are most interesting by Class specificity or expansion/loss within that Class.
However, you may also wonder which OGs contain gene families you are interested in a priori. 
Use FindContaining.py <path/to/OGs/fasta/dir> <file containing genes of interest, one per line>
Ranks OGs by the number of genes of interest they contain.

The genes trees provided by OrthoFinder are very useful for identifying Orthologs vs. Paralogs etc.
However, these can become confusing to the point of being unreadable when multiple species and genes are present.
To view the OG using a Network Graph, run 
$python OGfa2graph.py <input.fa> <output_directory>
Requires SpeClaD dictionary from before. Also requires a dictionary stating the unique header identifiers from each species assembly
e.g. Headers = {'ID1':'Species1', 'IDn':'Speciesn'}
Also requires Blast to be callable from the home directory. Make sure the module is loaded, or else edit the script. 
Runs a BLAST all v. all on each of the OG sequences and creates a pairwise.txt file of their relationships. 
Options for stringency of comparisons are available. 
This way, we can delve into when gene duplications occurred (phyogenetic relationships), which are Species/Class specific.
Colour coding the graph by species/class makes this much more comprehensible and the graph can be filtered according to similarity.nearest neighbour etc.

Finally, to obtain Chromosomal location of genes of interest for Synteny mapping, use:
$python Get_OG_member_loc.py <infile>
This requires the dictionary of unique identifiers used previously for each species. 
It outputs two resutls files for each species of interest: the gene names and their sequences specific to each species.
The list of gene names can then be run against a .gff annotaion file to give the chromosome, gene start, gene end and gene name.




