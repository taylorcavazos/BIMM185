### This script outputs all necessary files to fill genes and gene features tables.
### INPUTS: Genbank File and Organism for file naming
### OUTPUTS: A number of files used to fill SQL tables


import sys
from Bio import SeqIO # importing biopython
import gzip # importing library to read zipped files

# list with genomes you want to add to the database
gb_file = ["/home/linux/ieng6/bm185s/tcavazos/genomes/E_coli_K12_MG1655/GCF_000005845.2_ASM584v2_genomic.gbff.gz","/home/linux/ieng6/bm185s/tcavazos/genomes/A_tumafaciens/GCF_000576515.1_ASM57651v1/GCF_000576515.1_ASM57651v1_genomic.gbff.gz"]

# out files to write to different tables
out_genome = open("genomes.txt", "w")
out_replicon = open("replicons.txt", "w")
out_gene = open("genes.txt", "w")
out_exon = open("exons.txt", "w")
out_ext = open("extref.txt", "w")
out_syn = open("synonyms.txt", "w")
out_fun = open("functions.txt", "w")
# counters to maintain the current gene id, replicon id, and genome id for sql tables
# these fields should be updated to show the last id used for each field
genome_id, replicon_id, gene_id = 0,0,0

# begin looping through genomes
for i in range(0, len(gb_file)):
  genome_id = genome_id + 1 #increment genome id for every genome
  gnm_size, gnm_num_genes, num_reps = 0,0,0 # set counters for genome size, number of genes, and number of replicons to zero
  # open gzipped genome genbank file and parse through replicons in the genome
  with gzip.open(gb_file[i]) as f:
    for record in SeqIO.parse( f, "gb"):
      # getting genome and replicon information
      rep_shape = record.annotations.get("topology"); replicon_name = record.description; accession = record.name; date = record.annotations.get("date"); assembly = record.dbxrefs[2].split(":")[1]; rep_size = len(record.seq)
      # incrementing replicon counts and genome size
      num_reps = num_reps+1; replicon_id = replicon_id+1; num_genes_rep=0; gnm_size = gnm_size + rep_size
      # getting replicon type
      if "plasmid" in record.description.lower(): rep_type = "plasmid"
      else: rep_type = "chromosome"
      # extracting genome taxonomy
      domain = record.annotations["taxonomy"][0]
      if domain == "Bacteria": domain = "bacteria"
      elif domain == "Archaeon": domain = "archea"
      elif domain == "Eukaryota": domain = "eukarya"
      # looping through genes of replicon
      for feat in  record.features:
        #getting genome tax_id and name from source
	if feat.type=="source": 
          tax_id = feat.qualifiers.get("db_xref")[0].split(":")[1]
          genome_name = feat.qualifiers.get("organism")
        # only examine features if CDS, exclude pseudo genes 
	if feat.type=="CDS" and "pseudo" not in feat.qualifiers: 
	  # increment gene counters
          gene_id = gene_id+1; gnm_num_genes = gnm_num_genes+1; num_genes_rep = num_genes_rep+1 
	  # Extract protein id (external reference)
          if feat.qualifiers.get("protein_id") != None: 
	    protein_id = feat.qualifiers["protein_id"][0].split(".")[0]
	    out_ext.write(str(gene_id) + "\t" + "refseq" + "\t" + protein_id + "\n")
        # Extract external references
          if feat.qualifiers.get("db_xref")!=None:
	    for ref in feat.qualifiers.get("db_xref"):
              x_db = ref.split(":")[0]
              x_id = ref.split(":")[1]
	      out_ext.write(str(gene_id)+ "\t" + x_db + "\t" + x_id + "\n")
	# Extract gene synonyms
          if feat.qualifiers.get("gene_synonym") !=None:
	    for syn in feat.qualifiers.get("gene_synonym")[0].split("; "):
	      out_syn.write(str(gene_id) + "\t" + syn + "\n")
        # Extract locus tag
          locus_tag = feat.qualifiers["locus_tag"][0]
	# Extract gene function
          if feat.qualifiers.get("function") != None:
	    func = feat.qualifiers["function"][0]
	    out_fun.write(str(gene_id) + "\t" + func + "\n")
        # Extract gene name
	  if feat.qualifiers.get("gene") != None: name = feat.qualifiers["gene"][0]
	  elif feat.qualifiers.get("old_locus_tag") != None: name = feat.qualifiers["old_locus_tag"][0]
          else: name = " "
        # Extract direction of strand
          if feat.location.strand == 1: strand = "F"
          else: strand = "R"
        # Extract length of gene and number of exons
          num_exons = len(feat.location.parts)
          len_bp = len(feat)
        # Extract exon information such as start, stop, and length
          for l in feat.location.parts:
            left = l.start 
            right = l.end
            out_exon.write(str(gene_id) + "\t" + name + "\t" + str(left) + "\t" + str(right) + "\t" + str(len_bp) + "\n")
        # Extract product name of gene
          if feat.qualifiers.get("product") == None: product_name = " "
          else: product_name = feat.qualifiers["product"][0]
          # Write to genes table file 
          out_gene.write(str(gene_id) + "\t" + str(genome_id) + "\t" + str(replicon_id) + "\t" + str(locus_tag) +"\t"+protein_id+ "\t" + str(name) + "\t" + str(strand) + "\t" + str(num_exons) + "\t" + str(len_bp) + "\t"+ str(product_name) + "\n")
      # Write to replicons table file
      out_replicon.write(str(replicon_id) + "\t" + str(genome_id) + "\t" + str(replicon_name) + "\t" + str(rep_type) + "\t" + str(rep_shape) + "\t" + str(num_genes_rep) + "\t" + str(rep_size) + "\t" + str(accession) + "\t" + str(date) + "\n" )
  # Write to genomes table file 
  out_genome.write(str(genome_id) + "\t" + str(genome_name[0]) + "\t" + str(tax_id) + "\t" + domain + "\t" + str(num_reps) + "\t" + str(gnm_num_genes) + "\t" + str(gnm_size) + "\t" + str(assembly) + "\n")
