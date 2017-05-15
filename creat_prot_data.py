### This function takes an input of protein sequences and creates a protein data base that can be used for other downstream analysis
### The protein sequence number and identification number along with the sequence are outputted
 
import sys
import re

prot_file = sys.argv[1]
prot_data = open(prot_file).read().splitlines()

output_file = open('/home/linux/ieng6/bm185s/tcavazos/my_files/prot_database.txt', 'w')

prot_fixed = []
i = 0

while(i < len(prot_data)):
	if prot_data[i].startswith('>') == True:
		head = prot_data[i]
		i = i+1
	else:
		seq = ''
		while(i < len(prot_data) and prot_data[i].startswith('>') != True):
			seq = seq + prot_data[i]
			i = i+1
		prot_fixed.append((head, seq))

for i in range(0, len(prot_fixed)):
	match = re.search('>gnl\|[A-Z-]+\|([0-9]+)\|([1-9.A-Z]+)', prot_fixed[i][0])
	
	if match:
		output_file.write(str(match.group(2)+ "-"+ match.group(1)+ '\n'))
		output_file.write(str(prot_fixed[i][1]+"\n"))
	
	
