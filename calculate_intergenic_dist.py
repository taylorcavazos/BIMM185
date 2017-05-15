# Function that calculates intergenic distances between genes in an operon set
import sys
import MySQLdb
import re
# Connect to sql database
db = MySQLdb.connect()
# Create cursor to query tables
c = db.cursor()
# Open file for writing distances
out = open("intergenic_distances.txt", "w")
out.write("geneid_1\tgeneid_2\tdistance\n")

operons = open("operons.txt").read().splitlines() # list of operons and genes
products = open("gene_bnumber.txt").read() # list of genes and bnumbers

distances = [] # list to track distances between two adjacent genes
# Loop through all operons
for i in range(0, len(operons)):
  op = operons[i].split("\t")
  genes = op[1].split(",")
  # Only calculate distances if there is more than one gene in the operon
  if len(genes) > 1:
    # Loop through genes within operon
    for i in range(0, len(genes)-1):
      # search for bnumbers of two adjacent operons 
      s1 = re.search(str(genes[i])+"\t"+"(b[0-9]+)", products)
      s2 = re.search(str(genes[i+1])+"\t"+"(b[0-9]+)", products)
      # if there bnumbers were found save them
      if s1 and s2:
        bnum_1 = s1.group(1)
        bnum_2 = s2.group(1)
        # query genes table and extract gene id and strand info for both genes
        sql = "SELECT gene_id,strand FROM genes WHERE locus_tag='" + bnum_1 + "'"
        c.execute(sql)
        result_1 = c.fetchone()
        sql = "SELECT gene_id,strand FROM genes WHERE locus_tag='" + bnum_2 + "'"
        c.execute(sql)
        result_2 = c.fetchone()
        # if genes were found in the database extract exon info
        if result_1 != None and result_2 != None:
          id_1 = int(result_1[0]); strand_1 = result_1[1]
          id_2 = int(result_2[0]); strand_2 = result_2[1]
          # if the genes are in the forward strand calculate the distance by
          # left position gene 2 minus right position gene 1
          if strand_1 == "F" and strand_1 == strand_2:
            # get right position of gene 1
            sql = "SELECT right_pos FROM exons WHERE gene_id=" + str(id_1)
            c.execute(sql)
            right = int(c.fetchone()[0])
            # get left position of gene 2
            sql = "SELECT left_pos FROM exons WHERE gene_id=" + str(id_2)
            c.execute(sql)
            left = int(c.fetchone()[0])
            # write gene pair and distance to out file
            out.write(str(id_1) + "\t" + str(id_2) + "\t" + str(left - right) + "\n")
          # if the genes are in the reverse strand calculate the distance by
          # left position gene 1 minus right position gene 2
          elif strand_1 == "R" and strand_1 == strand_2:
            # get right position from gene 2
            sql = "SELECT right_pos FROM exons WHERE gene_id=" + str(id_2)
            c.execute(sql)
            right = int(c.fetchone()[0])
            # get left position from gene 1
            sql = "SELECT left_pos FROM exons WHERE gene_id=" + str(id_1)
            c.execute(sql)
            left = int(c.fetchone()[0])
            # write gene pair and distance to out file
            out.write(str(id_2) + "\t" + str(id_1) + "\t" + str(left-right) + "\n")
