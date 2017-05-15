### Function that uses biopython to parse a GenBank file to extract useful information
### INPUT: GenBank file
### OUTPUT: a tab seperated file conaining info such as accession, coordinates, strand, gene name, and more. 

from Bio import SeqIO # importing biopython
import gzip # importing library to read zipped files

# Open an output file for writing
out = open("GCF_000005845.2_ASM584v2_genomic_features.txt", "w")
# Writing the header to the file
out.write("Accession" + "\t"+"Coordinates"+ "\t"+ "Strand"+"\t"+"Gene Name"+ "\t"+"Locus Tag"+"\t"+"Synonyms"+ "\t"+"Protein Name"+"\t"+"Tax ID"+"\t"+"EC-number(s)"+"\t"+"External reference"+ "\n" )

# Opening the zipped GenBank file using gzip
with gzip.open("E_coli_K12_MG1655/GCF_000005845.2_ASM584v2_genomic.gbff.gz", "rt") as f:
  gb_record = SeqIO.read( f, "gb") # obtaining info from genbank file
  source = gb_record.features[0] # getting source tax ID from file
  tax_ID = source.qualifiers["db_xref"][0].split(":")[1]
  # Looping through features in genebank file 
  for feat in  gb_record.features:
    if feat.type == 'CDS': #keep features if they are CDS
       # output protein id to file if protein coding gene and psuedo if not for pseudogene
       if feat.qualifiers.get('protein_id') != None:
         out.write(feat.qualifiers['protein_id'][0]+ "\t")
       else:
         out.write("pseudo"+ "\n")
       # if there are multiple coordinates for the CDS, loop through them and separate by commas
       for l in feat.location.parts:
         out.write(str(l.start) + ":" + str(l.end) + ",")
       out.write("\t")
       # write what strand the gene is on
       out.write(str(feat.location.strand) + "\t")
       # write gene name for feature
       out.write(feat.qualifiers["gene"][0] + "\t")
       # write locus tag for feature
       out.write(feat.qualifiers["locus_tag"][0] + "\t")
       # write synonyms for feature
       out.write(feat.qualifiers["gene_synonym"][0] + "\t")
       
       # if there is product information, output it and if not output "-"
       if feat.qualifiers.get("product") != None:
         out.write(feat.qualifiers["product"][0] + "\t")
       else:
         out.write("-" + "\t")
       # write the tax ID obtained from the source
       out.write(tax_ID + "\t")
       # write all EC numbers
       if feat.qualifiers.get("EC_number") != None:
         for e in feat.qualifiers["EC_number"]:
            out.write(e + ",")
         out.write("\t")
       else:
         out.write("-" + "\t")
       # write the external reference to the file
       out.write(feat.qualifiers["db_xref"][0] + "\t")
    out.write("\n")
