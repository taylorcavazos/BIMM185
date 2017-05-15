# Function that outputs the distances between operons
import sys
import MySQLdb
import re
# Connect to MySQL database
db = MySQLdb.connect()
# Create cursor for database querying
c = db.cursor()

operons = open("operons.txt").read().splitlines() # list of operons and genes
products = open("gene_bnumber.txt").read() # list of genes and there bnumber

borders_fwd = {} # dictionary to contain forward strand genes on operon borders
borders_rev = {} # dictionary to contain reverse strand genes on operon borders
for ls in [operons]:
  # loop through operons
  for i in range(0, len(ls)):
    op = ls[i].split("\t") 
    # extract operon name and its genes
    name = op[0]; genes = op[1].split(",")
    # get genes on left and right border of operon
    gene_l = genes[0]; gene_r = genes[len(genes)-1]
    # search for bnumbers of border genes
    sl = re.search(str(gene_l)+"\t"+"(b[0-9]+)", products)
    sr = re.search(str(gene_r)+"\t"+"(b[0-9]+)", products)
    # if there bnumbers were found, find the gene id and strand 
    if sl and sr: 
      bnum_l = sl.group(1)
      bnum_r = sr.group(1)
      # query genes table for gene id and strand using bnumber
      sql = "SELECT gene_id,strand FROM genes WHERE locus_tag='" + bnum_l + "'"
      c.execute(sql)
      result_l = c.fetchone() 
      sql = "SELECT gene_id,strand FROM genes WHERE locus_tag='" + bnum_r + "'"
      c.execute(sql)
      result_r = c.fetchone()
      # if gene id was found for both genes extract the info
      if result_l != None and result_r != None:
        id_l = int(result_l[0]); strand_l = result_l[1]
        id_r = int(result_r[0]); strand_r = result_r[1]
        # Find left position of left border gene from exon table using gene id
        sql = "SELECT left_pos FROM exons WHERE gene_id=" + str(id_l)
        c.execute(sql)
        left_pos = int(c.fetchone()[0])
        # Find right position of right border gene from exon table using gene id
        sql = "SELECT right_pos FROM exons WHERE gene_id=" + str(id_r)
        c.execute(sql)
        right_pos = int(c.fetchone()[0])
        # If genes are on forward strand add them to forward dictionary
        if strand_l == "F" and strand_l == strand_r:
          if borders_fwd.get(name) == None: 
            borders_fwd[name] = (id_l, id_r, left_pos, right_pos)
        # If genes are on reverse strand add them to reverse dictionary
        elif strand_l == "R" and strand_l == strand_r:
          if borders_rev.get(name) == None:
            if len(genes) > 1: 
              borders_rev[name] = (id_r, id_l, right_pos,left_pos)
            else:
              borders_rev[name] = (id_l, id_r, left_pos, right_pos)
# Sort genes based on there left position
sorted_fwd = sorted(borders_fwd.items(), key=lambda x: x[1][2])
sorted_rev = sorted(borders_rev.items(), key=lambda x: x[1][2])
# Open output file to write distances between operons
distances = open("operon_distances.txt", "w")
# Loop through operon pairs in forward strand
for i in range(0, len(sorted_fwd)-1):
  op1_name = sorted_fwd[i][0]; op2_name = sorted_fwd[i+1][0]
  right_op1 = sorted_fwd[i][1][3]; left_op2 =  sorted_fwd[i+1][1][2]
  dist = left_op2-right_op1
  # write names of operon pair and distances between them to file
  distances.write(op1_name + "\t" + op2_name+ "\t" + "+" + "\t" + str(dist)+ "\n")
# Loop through operon pairs in reverse strand
for i in range(0, len(sorted_rev)-1):
  op1_name = sorted_rev[i][0]; op2_name = sorted_rev[i+1][0]
  right_op1 = sorted_rev[i][1][3]; left_op2 =  sorted_rev[i+1][1][2]
  dist = left_op2-right_op1
  # write names of operon pair and distances between them to file
  distances.write(op1_name + "\t" + op2_name +"\t" + "-" + "\t" + str(dist)+ "\n")
