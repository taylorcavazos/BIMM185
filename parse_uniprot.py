### Function that parses a UniProt file and outputs the tax id, organism, and taxonomy for every record
from Bio import SwissProt
import gzip
from collections import OrderedDict
# opening output file for writing
out = open("uniprot_sprot_archaea_features.txt", "w")
# maintain dict to hold lines so no duplicates are written to out file
lines_dict = OrderedDict()
#unzip uniprot file and use swissprot to parse file
with gzip.open("uniprot_sprot_archaea.dat.gz","rt") as f:
  for record in SwissProt.parse(f): # parsing file
    key = ",".join(record.taxonomy_id)
    value = str(record.organism) + "\t" + ";".join(record.organism_classification) + "\n"
    if lines_dict.get(key) == None:
      lines_dict[key] = value
# write uniq dict items to the file
out.write("NCBI Tax ID\tOrganism\tTaxomony\n")
for k, v in lines_dict.items():
  out.write(k + "\t" + v)
