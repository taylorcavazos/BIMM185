### This function calculates the CUI for each gene 
### INPUT: a file containing the bnumber, codon count, and length of gene
### OUTPUT: a file with bnumber and CUI for each gene

import sys

### Read in input file
codon_count = sys.argv[1]
codon_count = open(codon_count).read().splitlines()

# Open an output file
out = open("gene_CUI.txt", "w")

# Extract the total counts row from the file 
total_row = codon_count[len(codon_count)-1].split('\t')
total_codons_genome = int(total_row[len(total_row)-1])

### Calculate the CUI for each gene 
for i in range(1, len(codon_count)-1):
	line_split = codon_count[i].split('\t')
	gene_name = line_split[0]
	gene_length = int(line_split[len(line_split)-1])
	CUI = 0 #variable that will sum the CUI for each codon
	for j in range(1,len(line_split)-1):
		#q_c is the relative frequency of codon c in gene i
		q_c = float(line_split[j])/(gene_length/3) 
		#p_c is the probability of codon c in the genome
		p_c = float(total_row[j])/total_codons_genome
		CUI = CUI + float(q_c*p_c)
	out.write(gene_name + "\t" + str(CUI) + "\n")

