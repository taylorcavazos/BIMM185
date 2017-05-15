# Function that groups proteins that are related and sorts their scores in order from highest to lowest score
# Also a file with scores for each protein in ascending order is outputted

import sys
import bz2

prot_file = sys.argv[1]
prot_file_unzip = bz2.BZ2File(prot_file)
prot_file_read = prot_file_unzip.read().splitlines()

prot_dict = {}
count = 0
for i in range(0,len(prot_file_read)):
        curr_protein = prot_file_read[i].split("\t")
        if len(prot_dict) > 2000:
                break
        elif prot_dict.get(curr_protein[0]) == None:
                prot_dict[curr_protein[0]] = [(curr_protein[1],curr_protein[3])]
        else:
                prot_dict[curr_protein[0]].append((curr_protein[1],curr_protein[3]))

### Write dictionary to a file
output_file = open('/home/linux/ieng6/bm185s/tcavazos/my_files/prot_relations.txt', 'w')
for k, v in prot_dict.items():
        for val in sorted(v, key=lambda x:x[1], reverse=True) :
                output_file.write(str(k+'\t'+val[0]+'\t'+val[1]+'\n'))

#### find protein with max number of related proteins

max_score = float("-inf")
prot_max = None
for keys, values in prot_dict.items():
           if len(values) > max_score:
                max_score = len(values)
                prot_max = keys
print(prot_max)
print(max_score)
