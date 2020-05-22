#!/usr/bin/env python

#Usage: python OGfa2graph.py <input.fa> <output_directory>
#e.g. python OGfa2graph.py OGs_Results_May16/Bivalvia_Dups/OG0000012.fa OG0000012_HSP -p 40 -l 50
"""
Takes an Orthogroup.fasta file and runs Blastp all versus all to generate a pairwise .txt file.
This can be used to generate a graph using Graphia. Nodes are also annotated according to species and class. 
Requires Blastp module to be loaded, and the ability to create a megablastDB
e.g. module load igmm/apps/blast/2.2.26
"""
import os, glob, commands
import string
import sys, optparse
import math
import shlex
import subprocess
from operator import itemgetter

Headers = {'MDBACRIN':'Acanthochitona_crinita', 'Afu':'Achatina_fulica', \
'A_californica_':'Aplysia_californica', "Bpl_scaf":'Bathymodiolus_platifrons', 'B_glabrata_':'Biomphalaria_glabrata', \
"C_gigas_":'Crassostrea_gigas', "C_virginica_":'Crassostrea_virginica', "MDBCPLIC":'Cristaria_plicata',\
 'RUS':'Elysia_chlorotica', 'MDBGTOLM':'Gadila_tolmiei', 'MDBGPELL':'Gymnomenia_pellucida', \
 'MDBLHYAL':'Laevipilina_hyalina', 'Laternula_elliptica_':"Laternula_elliptica", 'Lotgi':'Lottia_gigantea', 'L_gigantea_':'Lottia_gigantea',\
 'L_stagnalis_':'Lymnaea_stagnalis', 'M_yessoensis_':'Mizuhopecten_yessoensis', "Mph_scaf":'Modiolus_philippinarum',\
 'MDBMAREN':'Mya_arenaria', 'Mya_truncata_':'Mya_truncata', 'MDBMEDUL':'Mytilus_edulis', \
 'MDBMGALL':'Mytilus_galloprovincialis', 'od_comp':'Octopoteuthis_deletron', 'O_bimaculoides_':'Octopus_bimaculoides', \
 "O_vulgaris_":'Octopus_vulgaris', "MDBPMAXI":'Pecten_maximus', \
 'pfu_':'Pinctada_fucata', "P_canaliculata_":'Pomacea_canaliculata', 'sv_':'Scutopus_ventrolineatus', \
 'vi_':'Vampyroteuthis_infernalis', 'MDBWARGE':'Wirenia_argentea'}
#This dict yields the species associated with a given pattern in a Fasta header - Pattern:Species

TaxClass_Colours = {"Bivalvia":"#1f77b4", "Caudofoveata":"#ff7f0e", "Cephalapoda":"#2ca02c", "Gastropoda":"#d62728", \
"Monoplacophora":"#9467bd", "Polyplacophora":"#8c564b", "Scaphopoda":"#e377c2", "Solenogastres":"#7f7f7f"}

SpeClaD = {'Acanthochitona_crinita': 'Polyplacophora', 'Achatina_fulica': 'Gastropoda', 'Aplysia_californica': 'Gastropoda', \
'Bathymodiolus_platifrons': 'Bivalvia', 'Biomphalaria_glabrata': 'Gastropoda', 'Crassostrea_gigas': 'Bivalvia', \
'Crassostrea_virginica': 'Bivalvia', 'Cristaria_plicata': 'Bivalvia', 'Elysia_chlorotica': 'Gastropoda', \
'Gadila_tolmiei': 'Scaphopoda', 'Gymnomenia_pellucida': 'Solenogastres', 'Laevipilina_hyalina': 'Monoplacophora', \
'Laternula_elliptica': 'Bivalvia', 'Lottia_gigantea': 'Gastropoda', 'Lymnaea_stagnalis': 'Gastropoda', \
'Mizuhopecten_yessoensis': 'Bivalvia', 'Modiolus_philippinarum': 'Bivalvia', 'Mya_arenaria': 'Bivalvia', \
'Mya_truncata': 'Bivalvia', 'Mytilus_edulis': 'Bivalvia', 'Mytilus_galloprovincialis': 'Bivalvia', \
'Octopoteuthis_deletron': 'Cephalapoda', 'Octopus_bimaculoides': 'Cephalapoda', 'Octopus_sinensis': 'Cephalapoda', \
'Octopus_vulgaris': 'Cephalapoda', 'Pecten_maximus': 'Bivalvia', 'Pinctada_fucata': 'Bivalvia', \
'Pomacea_canaliculata': 'Gastropoda', 'Scutopus_ventrolineatus': 'Caudofoveata', 'Vampyroteuthis_infernalis': 'Cephalapoda', \
'Wirenia_argentea': 'Solenogastres'}

def megablast2ncol(fastaFile, alignFile, outFile, matchaccuracy, matchlength):
  """
  Transform megablast output in ncol format and filter out hits with a match accuracy and a matchlength below the given threshold
  The script reads the output of megablast (-D 3 option for output format).
  It takes a lot of memory if the megablast output file is big.
  """
  fastaF = open(fastaFile, "r")
  alignF = open(alignFile, "r")
  outF = open(outFile, "w")

  dicFasta = {}
  seqA = set()
  for line in fastaF:
    if(line[0] == ">"):
      seqname = string.split(string.split(line, ">")[1],"\n")[0]
      seqname = seqname.split()[0]
      seqA.add(seqname)
    else:
      seq = string.split(line,"\n")[0]
      if seqname in dicFasta:
        dicFasta[seqname] += len(seq)
      else:
        dicFasta[seqname] = len(seq)
  fastaF.close()

  #read output of megablast (-D 3 option for output format has been used)
  dicEdges = {}
  edgesList = []

  for line in alignF:
    if (line[0] != "#"):
      l = line.split()
      q = l[0]
      t = l[1]
      if (t != q) and (t != 'Query_1') and (q != 'Query_1'): #Query_1 added after weird megablast error...this may omit some relationships but can't figure it out
        iden = float(l[2])/100
        align_len = int(l[3])
        mismatches = int(l[4])
        bitScore = int(float(l[11])) #round off to int
        q_len = dicFasta[q]
        t_len = dicFasta[t]
        
        coverage = align_len/float(q_len)
        
        if (iden >= matchaccuracy) and (coverage >= matchlength):
          try:
            dicEdges[q,t]
            dicEdges[t,q]
          except:
            dicEdges[q,t] = bitScore
            dicEdges[t,q] = bitScore
            edgesList.append([q,t,bitScore])
  alignF.close()
  del(dicFasta)
  del(dicEdges)

  edgesList_sorted = sorted(edgesList, key=itemgetter(2), reverse=True)
  del(edgesList)
  
  totNode = set()
  for edge in edgesList_sorted:
    outF.write(edge[0] + "\t" + edge[1] + "\t" + str(edge[2]) + "\n")
    totNode.add(edge[0])
    totNode.add(edge[1])
  outF.close()

  del(edgesList_sorted)

  print len(totNode)


if __name__ == "__main__":

  scriptDir = os.path.dirname(os.path.realpath(__file__))
  print "Script directory is " + scriptDir

  parser = optparse.OptionParser("usage: %prog [options] fastaFile outputDir")
  parser.add_option("-p", "--percentage", dest="similarity", default=75, type="int", help="specify sequence similarity percentage. Default=85.")
  parser.add_option("-W", "--wordsize", dest="wordsize", default=2, type="int", help="specify megablast word size. Default=6.")
  parser.add_option("-a", "--parallel", dest="numthreads", default=8, type="int", help="specify megablast number of threads. Default=1.")
  parser.add_option("-l", "--coverage", dest="coverage", default=55, type="int", help="specify percentage of length coverage for creating the graph. Default=55.")
  parser.add_option("-c", "--contigassembly", action="store_true", dest="contig", default=False, help="create contings using CAP3 and create BioLayout class files. Default=False.")
  parser.add_option("-k", "--keepoutput", action="store_true", dest="keepoutput", default=False, help="keep intermediate data files. Default=False.")
  parser.add_option("-q", "--verbose", action="store_false", dest="verbose", default=True, help="quit verbose. Default=True.")
  (options, args) = parser.parse_args()

  if (len(args) != 2):
    parser.error("incorrect number of arguments")

  fastaFile = os.path.abspath(args[0])
  outputDir = os.path.abspath(args[1])
  percentage = options.similarity
  wordsize = options.wordsize
  numthreads = options.numthreads
  coverage = options.coverage
  contig = options.contig
  keepoutput = options.keepoutput
  verbose=options.verbose

  if not os.path.isdir(outputDir):
    os.makedirs(outputDir)
  else:
    print "Output directory already exists!"
    sys.exit(1)

  os.chdir(outputDir)
  nameDB = os.path.splitext(os.path.basename(fastaFile))[0]
  cmd = "makeblastdb -dbtype prot -out " + nameDB + " -input_type fasta -in " + fastaFile + " -title " + nameDB + " -max_file_sz 2GB"
  if verbose:
    print cmd
  args = shlex.split(cmd)
  p = subprocess.Popen(args, stdout=subprocess.PIPE)

  output = p.communicate()[0]
  if p.returncode != 0:
    print output
    print "...failed"
    exit(p.returncode)

  numLines = sum(1 for line in open(fastaFile))
  linesPerThread = -(-numLines // numthreads)

  # We need an even number of lines because fasta format is pairs of lines
  if linesPerThread % 2 != 0:
    linesPerThread += 1

  remainingLines = 0
  thread = 0
  fastaThreadFile = None
  for line in open(fastaFile):
    fastaThreadFilename = outputDir + "/" + nameDB + ".t" + `int(thread)` + ".fasta"
    if remainingLines == 0:
      if fastaThreadFile != None:
        fastaThreadFile.close()

      fastaThreadFile = open(fastaThreadFilename, 'a+')
      remainingLines = linesPerThread
      thread += 1

    fastaThreadFile.write(line)
    remainingLines -= 1

  if fastaThreadFile != None:
    fastaThreadFile.close()

  # May not be able to divide into requested number of threads
  numthreads = thread

  megablastProcesses = []
  for thread in range(0, numthreads):
    fastaThreadFilename = outputDir + "/" + nameDB + ".t" + `int(thread)` + ".fasta"
    megablastThreadFilename = outputDir + "/" + nameDB + ".t" + `int(thread)` + "_megablast.txt"

    cmd = "blastp -db " + nameDB + " -query " + fastaThreadFilename + " " \
        "-qcov_hsp_perc " + str(percentage) + " -word_size " + str(wordsize) + " " \
        "-xdrop_gap 40 -num_alignments 90000000 -soft_masking true -lcase_masking " + \
        "-outfmt \"7 qacc sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" " + \
        "-out " + megablastThreadFilename
    if verbose:
      print cmd

    args = shlex.split(cmd)
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    megablastProcesses.append({'process' : p, 'blastFilename' : megablastThreadFilename, 'fastaFilename' : fastaThreadFilename})

  megablastFilename = outputDir + "/" + nameDB + "_megablast.txt"
  with open(megablastFilename, 'w') as megablastFile:
    for megablastProcess in megablastProcesses:
      p = megablastProcess['process']
      megablastThreadFilename = megablastProcess['blastFilename']
      fastaThreadFilename = megablastProcess['fastaFilename']

      output = p.communicate()[0]
      if p.returncode != 0:
        print output
        print "...failed"
        exit(p.returncode)

      with open(megablastThreadFilename) as megablastThreadFile:
        for line in megablastThreadFile:
          megablastFile.write(line)
      if not keepoutput:
        os.remove(megablastThreadFilename)
        os.remove(fastaThreadFilename)

  pairwise_txtFile = outputDir + "/" + nameDB + "_pairwise.txt"
  

  matchaccuracy = percentage/float(100)
  matchlength = coverage/float(100)
  megablast2ncol(fastaFile, megablastFilename, pairwise_txtFile, matchaccuracy, matchlength)

  execCmd = os.popen(cmd, "r")
  outCmd = execCmd.read()
  execCmd.close()

  if not keepoutput:
    os.remove(megablastFilename)
#Add Class/Species annotation to the bottom of pairwise_txtFile
genes = []
with open(fastaFile,'r') as OGfa:
    for line in OGfa:
        if(line[0] == ">"):
            seqname = (line.strip(">")).strip()
            genes = genes + [seqname]

with open(pairwise_txtFile, 'a') as ClassAssign:
  ClassAssign.write("//NODECLASS\t//Gene\t//ID\t//Grouping\n")
  for gene in genes:
    species = "".join([spec for pattern, spec in Headers.items() if pattern in gene])
    Class = SpeClaD[species]
    lines = "//NODECLASS\t//%s\t//%s\t//Species\n//NODECLASS\t//%s\t//%s\t//Class\n" % (gene, species, gene, Class)
    ClassAssign.write(lines)
  for key in TaxClass_Colours:
    colour_line = "//NODECLASSCOLOR\t//%s\t//Class\t//%s\n" % (key, TaxClass_Colours[key])
    ClassAssign.write(colour_line)
