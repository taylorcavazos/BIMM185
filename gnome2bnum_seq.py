### Function that takes in a genome annotation in tabular format 
### and the genome sequences in fasta format
### INPUTS: annotation table and whole genome
### OUTPUTS: a file containing the b number and DNA sequence seperated by a tab
 
import sys 
import csv
import textwrap

DNA_dict = {'A': 'T', 'C': 'G', 'T':'A', 'G':'C'}
file_out = open('fasta_2.txt', 'w')

### read in the table
table = sys.argv[1]
with open(table) as table:
                reader = csv.reader(table, delimiter="\t")
                table_list = list(reader)

### read in the genome sequence
genome = sys.argv[2] 
genome = open(genome,'r').read().splitlines()
genome_str = ''.join(genome[1:])

### loop through genes
for i in range(1, len(table_list)-1):
	gene_id = table_list[i][7]
	start = int(table_list[i][2])
	stop = int(table_list[i][3])
	strand = table_list[i][4]
	file_out.write(str(gene_id+ '\t'))
	seq = genome_str[start-1:stop]
	if strand == "+":
		file_out.write(seq)
		file_out.write("\n")
	### if on reverse strand, find the reverse complement
	else:                                         
		rev_seq = seq[::-1]
		rev_comp = ''
		for r in rev_seq:
			rev_comp = rev_comp + DNA_dict.get(r)	
		file_out.write(rev_comp)
                file_out.write('\n')			
