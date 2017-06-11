# Pipeline that uses a GEO dataset and the output of GSEA to analyze the results

# Import necessary libraries for analysis
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import decimal
import MySQLdb

# Read in GEO dataset GSE43837
series = open("GSE43837_series_matrix.txt").read().splitlines()

# Step 1: Prepare input data for GSEA
# Find position in file where series matrix starts, ignore all data above this position
rows = []
pos = 0
for line in series:
  if re.search("!Sample_source_name_ch1", line):
    phenotypes = line.split("!Sample_source_name_ch1\t")[1].split("\t")
  if re.search("!series_matrix_table_begin", line):
    break
  pos = pos+1
for i in range(pos+1, len(series)-1):
  rows.append(series[i])
# extract the number of samples from dataset
num_samples = len(rows[0].split("\t")[1:])
# output phenotype.cls file (input for GSEA)
phen = open("phenotype_GSE43837.cls", "w")
phen.write(str(num_samples) + "\t2\t1\n")
phen.write("# MET PRIMARY\n")
phen.write("\t".join(list(np.repeat(["0", "1"], [19, 19])))) 
phen.close()
# output gene_expression.txt file (input for GSEA)
expr = open("gene_expression_GSE43837.txt", "w")
header = rows[0].split("\t")
expr.write("NAME\tDESCRIPTION\t")
for i in range(1, len(header)-1):
  expr.write(str(eval(header[i])) + "\t")
expr.write(str(eval(header[len(header)-1])) + "\n")
for i in range(1, len(rows)):
  r = rows[i].split("\t")
  expr.write(str(eval(r[0]))+"\tna\t" + "\t".join(r[1:]) + "\n")
expr.close()

# Assuming GSEA was performed, gene set was chosen as signature, and core enriched genes in pathway were converted to their probes and included as an input file...
# Step 2: Analyze gene set in comparison to known signatures for HER2+ brain metastases samples
# Read in gene list and gene expression file
GSE43837 = pd.read_csv("gene_expression_GSE43837.txt", sep="\t", index_col=0)
DEK_probes = open("DEK_genes_probes_hgu133.txt").read().splitlines()
# connect to mySQL database
db = MySQLdb.connect()
# create cursor for database querying
c = db.cursor()
# remove unecessary description column
del GSE43837["DESCRIPTION"]
# make dictionary for probes and gene symbol
DEK = {}
for i in range(1, len(DEK_probes)):
    probe_gene = DEK_probes[i].split(" ")
    DEK[probe_gene[0][1:-1]] = probe_gene[1][1:-1]
# subset dataset to only include core enriched probes in DEK pathway
df = pd.DataFrame(columns=GSE43837.columns)
for k in DEK.keys():
    if k in list(GSE43837.index):
        df.loc[GSE43837.ix[k,:].name] = GSE43837.ix[k,:]

# plot the DEK oncogene expression signature for primary and metastasis samples
y = []
x = [1,2]
y_met = ()
y_pri = ()
for col in list(df.columns):
  # access mysql table consisting of patient, tumor type, HER2 status, ER status, and Age
  sql = "SELECT Tumor_Type FROM breast_patient_info WHERE Sample_ID='" + col + "'"
  c.execute(sql)
  tumor_type = c.fetchone()[0]
  if tumor_type == "brain metastasis":
    y_met = y_met + (df[[col]].mean().values[0],)
  if tumor_type == "primary breast tumor":
    y_pri = y_pri + (df[[col]].mean().values[0],)
y.append(y_met)
y.append(y_pri)
for xe, ye in zip(x, y):
  plt.scatter([xe]*len(ye), ye)
plt.xticks([1,2])
plt.boxplot([y_met, y_pri])
plt.axes().set_xticklabels(['Metastases', 'Primary'])
wilcox = round(decimal.Decimal(stats.wilcoxon(y_met, y_pri)[1]), 4)
plt.title("Wilcoxon p=" + str(wilcox),size=14)
plt.ylabel("DEK Oncogene Encriched Signature")
plt.show() 

# Plot the BRCA1 (probe = g2218153_3p_a_at) expression for samples
BRCA_exp = GSE43837.ix["g2218153_3p_a_at",:]
y = []
x = [1,2]
y_met = ()
y_pri = ()
for col in list(df.columns):
  sql = "SELECT Tumor_Type FROM breast_patient_info WHERE Sample_ID='" + col + "'"
  c.execute(sql)
  tumor_type = c.fetchone()[0]
  if tumor_type == "brain metastasis":
    y_met = y_met + (df[[col]].mean().values[0],)
  if tumor_type == "primary breast tumor":
    y_pri = y_pri + (df[[col]].mean().values[0],)
y.append(y_met)
y.append(y_pri)
for xe, ye in zip(x, y):
  plt.scatter([xe]*len(ye), ye)
plt.xticks([1,2])
plt.boxplot([y_met, y_pri])
plt.axes().set_xticklabels(['Metastases', 'Primary'])
wilcox = round(decimal.Decimal(stats.wilcoxon(y_met, y_pri)[1]), 4)
plt.title("Wilcoxon p=" + str(wilcox),size=14)
plt.ylabel("BRCA1 Expression (g2218153_3p_a_at)")
plt.show()

# Plot the BRCA1 (probe = g6552300_3p_a_at) expression for samples
BRCA_exp = GSE43837.ix["g6552300_3p_a_at",:]
y = []
x = [1,2]
y_met = ()
y_pri = ()
for col in list(df.columns):
  sql = "SELECT Tumor_Type FROM breast_patient_info WHERE Sample_ID='" + col + "'"
  c.execute(sql)
  tumor_type = c.fetchone()[0]
  if tumor_type == "brain metastasis":
    y_met = y_met + (df[[col]].mean().values[0],)
  if tumor_type == "primary breast tumor":
    y_pri = y_pri + (df[[col]].mean().values[0],)
y.append(y_met)
y.append(y_pri)
for xe, ye in zip(x, y):
  plt.scatter([xe]*len(ye), ye)
plt.xticks([1,2])
plt.boxplot([y_met, y_pri])
plt.axes().set_xticklabels(['Metastases', 'Primary'])
wilcox = round(decimal.Decimal(stats.wilcoxon(y_met, y_pri)[1]), 4)
plt.title("Wilcoxon p=" + str(wilcox),size=14)
plt.ylabel("BRCA1 Expression (g6552300_3p_a_at)")
plt.show()

# Explore the effect of ER status on DEK expression
ER_pos_met = []
ER_pos_pri = []
ER_neg_met = []
ER_neg_pri = []
for col in df.columns:
  sql = "SELECT Sample_ID, Tumor_Type, ER_status FROM breast_patient_info WHERE Sample_ID='" + col + "'"
  c.execute(sql)
  status = c.fetchone()
  if status[1] == "brain metastasis":   
    if status[2] == "+": 
      ER_pos_met.append(df[[status[0]]].mean().values[0])
    elif status[2] == "-":
      ER_neg_met.append(df[[status[0]]].mean().values[0])
  elif status[1] == "primary breast tumor":
    if status[2] == "+":
      ER_pos_pri.append(df[[status[0]]].mean().values[0])
    elif status[2] == "-":
      ER_neg_pri.append(df[[status[0]]].mean().values[0])
y = []
x = [1,2,3,4]
y.append(ER_neg_met)
y.append(ER_pos_pri)
y.append(ER_neg_pri)
for xe, ye in zip(x, y):
  plt.scatter([xe]*len(ye), ye)
plt.xticks([1,2,3,4])
plt.boxplot([ER_pos_met, ER_neg_met, ER_pos_pri, ER_neg_pri])
plt.axes().set_xticklabels(['Met ER+','Met ER-', 'Primary ER+', 'Primary ER-'])
plt.title("ER Status Versus DEK Signature Expression",size=14)
plt.ylabel("DEK Expression")
plt.show()
