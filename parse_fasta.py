### Function that parses a fasta file 
### INPUT: Fasta file
### OUTPUT: Accession and sequence of every protein

from Bio import SeqIO # importing SeqIO from biopython to parse file
import gzip # importing gzip to read zipped file

# opening output file to write to and writing the header
out = open("GCF_000005845.2_ASM584v2_protein_features", "w")
out.write("Accession"+ "\t"+"Protein Sequence"+ "\n")

with gzip.open("E_coli_K12_MG1655/GCF_000005845.2_ASM584v2_protein.faa.gz","rt") as f:
  fasta = SeqIO.parse(f, "fasta") # parsing file
  for prot in fasta: # looking through every line that was parsed
    accession, sequence = prot.id, str(prot.seq) # save accession and sequence for every protein using id and seq
    out.write(accession + "\t"+ sequence + "\n") # writing output to file
