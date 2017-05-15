### Function that estimates the frequencies of codons 
### INPUT: file containing bnumber gene separated by sequence
### OUTPUT: table of gene name codon frequency and length of gene

import sys
import textwrap
from collections import OrderedDict

### read in bnumber sequence file
gene_file = sys.argv[1]
genes_seq = open(gene_file).read().splitlines()

### open output file
output = open('codon_count_table.txt', 'w')

### list of 64 codons
codons = ['ATG', 'ATT', 'ATC', 'ATA', 'CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG',
'GTT', 'GTC', 'GTA', 'GTG', 'TTT', 'TTC', 'TGT', 'TGC', 'GCT', 'GCC', 'GCA', 'GCG',
'GGT', 'GGC', 'GGA', 'GGG', 'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 'ACC', 'ACA', 'ACG', 'TCT', 'TCC', 'TCA',
 'TCG', 'AGT', 'AGC', 'TAT', 'TAC', 'TGG', 'CAA', 'CAG', 'AAT', 'AAC', 'CAT', 'CAC', 'GAA', 'GAG',
'GAT', 'GAC', 'AAA', 'AAG', 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG', 'TAA', 'TAG', 'TGA']

### dictionary to make global counts
codon_global_count = OrderedDict()

### write header of table
output.write('Gene'+ '\t')
for c in codons:
	output.write(c + '\t')
	codon_global_count[c] = 0
output.write('Length'+ '\n')

### nested dictionary to keep track of individual gene codon counts
gene_codon_count = OrderedDict()
total_gene_len = 0
for genes in genes_seq:
	gene = genes.split('\t')
	length_gene = len(gene[1])
	total_gene_len = total_gene_len + length_gene
	### check to ensure genes are divisible by three 
	if length_gene % 3 == 0:
		gene_codon_count[gene[0]] = OrderedDict()
		for c in codons: 
			gene_codon_count[gene[0]][c]=0
		### split the gene into its codons
		seq_split = textwrap.wrap(gene[1], 3)
		for cod in seq_split:
			codon_global_count[cod] = codon_global_count.get(cod) + 1
			gene_codon_count[gene[0]][cod] = gene_codon_count[gene[0]][cod]+1
		### write the gene and the counts of each codon to the file
		output.write(str(gene[0] + '\t'))
		for k in gene_codon_count[gene[0]].keys():
			output.write(str(gene_codon_count[gene[0]][k])+ '\t')
		output.write(str(length_gene)+ '\n')

### write the total counts for each codon to file 
output.write('Totals' + '\t')
for k in codon_global_count.keys():
	output.write(str(codon_global_count.get(k)) + '\t')
output.write(str(total_gene_len))
output.write('\n')
