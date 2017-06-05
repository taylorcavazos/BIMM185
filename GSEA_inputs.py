# Function to convert matrix series files from NCBI GEO into GSE inputs
# INPUT: NCBI GEO series matrix
# OUTPUTS: (1) gene expression file (.txt) containing probes versus samples and expression values
########## (2) phenotype (.cls) file
# specifications for the above files can be found at (http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) 

import sys
import re 

# read in series data
series = open(sys.argv[1]).read().splitlines() #input series matrix from GEO and read 
rows = []
# find position in file where series matrix starts, ignore all data above this position
pos = 0
for line in series:
  if re.search("!Sample_source_name_ch1", line):
    phenotypes = line.split("!Sample_source_name_ch1\t")[1].split("\t")
  if re.search("!series_matrix_table_begin", line):
    break
  pos = pos+1    
for i in range(pos+1, len(series)-1):
  rows.append(series[i])
# extract the number of samples from dataset
num_samples = len(rows[0].split("\t")[1:])

# output phenotype .cls file
phen_file = open("phenotypes.cls", "w")
phen_file.write(str(num_samples) + "\t2\t1\n")
phen_file.write("# MET PRIMARY\n")
# write phenotype for each sample (0 for metastatic 1 for primary)
for p in range(0, len(phenotypes)-1):
  if phenotypes[p] == '"brain metastasis"':
    phen_file.write("0 ")
  elif phenotypes[p] == '"primary breast tumor"':
    phen_file.write("1 ")
if phenotypes[p+1] == '"brain metastasis"': phen_file.write("0\n")
else: phen_file.write("1\n")

# output expression data
expr = open("gene_expression.txt", "w") 
header = rows[0].split("\t")
expr.write("NAME\tDESCRIPTION\t")
for i in range(1, len(header)-1):
  expr.write(str(eval(header[i])) + "\t")
expr.write(str(eval(header[len(header)-1])) + "\n")
for i in range(1, len(rows)):
  r = rows[i].split("\t")
  expr.write(str(eval(r[0])) + "\tna\t" + "\t".join(r[1:]) + "\n") 


