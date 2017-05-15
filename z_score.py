###This function outputs the highest z-score for each report.tbl file in the sub directorys
### Input: file with directories pointing to report files for sub directories in dirtree

import sys
import csv

files = sys.argv[1]

output_file = open("output.txt", 'w')
files = open(files).read().splitlines()
for f in files:
	with open(f) as f_open:
    		reader = csv.reader(f_open, delimiter="\t")
    		d = list(reader)		
		
		output_file.write(f.split('/')[0]) ### write the directory name to the file
		output_file.write('\t')
		output_file.write(d[2][3]) ### write the highest z score to the file
		output_file.write('\n')





