from __future__ import division
import sys
import argparse
import time
start_time = time.time()

parser=argparse.ArgumentParser()
parser.add_argument("-i","--INPUT",help="The input file containing the stats for the gene types you wish to compare to the background.",type=str,required=True)
parser.add_argument("-o","--OUTPUT",help="The output generated containing resampled differences for each stat.",type=str,required=True)
parser.add_argument("-b","--BOOTSTRAP",help="Number of resamples to do for each category and each statistic.",type=int,default=1000)
parser.add_argument("-w","--WINDOWSIZE",help="The window size to sample genes near the focal gene across.",type=int, default=100000)
parser.add_argument("-s","--SIZEDIFFERENCE",help="The maximum difference in size (bp) between the focal gene and the sampled genes.",type = int, default = 10000)
parser.add_argument("-v","--VALUE",help="The stat or stats (comma seperated, or line by line in a text file) to resample, given as a column name.",required=True)
parser.add_argument("-c","--CLASS",help="The class or classes (comma seperated, or line by line in a text file) to find difference from the background.",required=True)
parser.add_argument('-quiet','--QUIET', help="If set, does not print any messages.", type=bool,default=False)
args=parser.parse_args()

"""
input should be a tab delimited bed file with a header in the following format:
chr	start	stop	class	group	STAT

Any number of columns with different statistics for each gene can be included after the class function.
Columns must be named as above, gene name is not required.
If a stat is going to be resampled, its column name should be called with -v
If a group column is present and multiple groups (e.g. species, populations) are found within the column, the script will iterate over each group seperately.
"""
import os
import pandas as pd
import numpy as np
from collections import defaultdict

input_file = pd.read_table(args.INPUT)

"check if column distinguishing classes of interest is included."
if 'class' not in input_file.columns:
	print("No class detected in input, exiting.")
	quit()

"check if all columns of statistics are present in dataframe. Ignore absent statistics."
values = []
append = values.append
if os.path.isfile(args.VALUE) == True:
	for line in open(args.VALUE):
		if line.rstrip() in input_file.columns:
			append(line.rstrip())
		else:
			print(i + " not a statistic in dataframe, ignoring.")
else:
	for i in args.VALUE.split(","):
		if i in input_file.columns:
			append(i)
		else:
			print(i + " not a statistic in dataframe, ignoring.")

mean = np.mean

windows = args.WINDOWSIZE
difference = args.SIZEDIFFERENCE

input_file['gene_size'] = input_file['end'] - input_file['start']

def filter_sample(rows,controls):
	controls_1 = controls[controls['chr'] == rows['chr']]
	controls_2 = controls_1[controls_1['start'] > rows['start'] - windows]
	controls_3 = controls_2[controls_2['start'] < rows['start'] + windows]
	controls_4 = controls_3[rows['gene_size'] * difference > controls_3['gene_size']]
	controls_5 = controls_4[controls_4['gene_size'] * difference > rows['gene_size']]
	return(controls_5)

if 'group' in input_file.columns:
	"Generate output table."
	output_file = pd.DataFrame(columns=["group","class","rep"])
	for i in values:
		output_file[i] = []
	append = output_file.append
	"Iterate over classes to resample differences and generate output."
	if len(set(list(input_file['group'].unique()))) > 1:
		if args.QUIET == False:
			print("Resampling multiple groups.")
	for xx in set(list(input_file['group'].unique())):
		input_file2 = input_file[input_file['group'] == xx]
		classes = []
		if os.path.isfile(args.CLASS) == True:
			for line in open(args.CLASS):
				if line.rstrip() in list(input_file2["class"].unique()):
					classes.append(line.rstrip())
				else:
					print(i + " not a class in " + xx + ", ignoring.")
		else:
			for i in args.CLASS.split(","):
				if i in list(input_file2["class"].unique()):
					classes.append(i)
				else:
					print(i + " not a class in " + xx + ", ignoring.")
		control_samples_1 = input_file2[~input_file2['class'].isin(classes)]
		for i in classes:
			case_samples = input_file2[input_file2['class'] == i]
			control_samples_2 = control_samples_1[control_samples_1['chr'].isin(list(case_samples['chr'].unique()))]
			count = 0
			iters = case_samples.iterrows
			for j in range(0,args.BOOTSTRAP):
				count += 1
				if args.QUIET == False:
					if count % (args.BOOTSTRAP/10) == 0:
						print("%s %s resamples processed in %s" % (str(j + 1),str(i),xx))
				stats = defaultdict(list)
				for index,row in iters():
					control_subset = filter_sample(row,control_samples_2)
					if len(control_subset) > 0:
						row = row[values]
						control_subset_2 = control_subset[values].sample()
						for k, v in dict(row-control_subset_2).items():
							stats[k].extend(v.tolist())
				vv = {"group":xx,"class":i,"rep":j}
				stats = dict(stats)
				for v in stats:
					vv[v] = mean(stats[v])
				df  = pd.DataFrame([vv], columns=vv.keys())
				output_file = pd.concat([output_file, df], axis =0,sort=False,ignore_index=True)

else:
	"Generate output table."
	output_file = pd.DataFrame(columns=["class","rep"])
	for i in values:
		output_file[i] = []
	append = output_file.append
	"check if all classes are present in dataframe. Ignore absent classes"
	classes = []
	append = classes.append
	if os.path.isfile(args.CLASS) == True:
		for line in open(args.CLASS):
			if line.rstrip() in list(input_file["class"].unique()):
				append(line.rstrip())
			else:
				print(i + " not a class in dataframe, ignoring.")
	else:
		for i in args.CLASS.split(","):
			if i in list(input_file["class"].unique()):
				append(i)
			else:
				print(i + " not a class in dataframe, ignoring.")
	"If no statistics or classes are chosen, exit."
	if len(classes) == 0 or len(values) == 0:
		print("No classes or values selected, exiting.")
		quit()
	control_samples_1 = input_file[~input_file['class'].isin(classes)]
	for i in classes:
		case_samples = input_file[input_file['class'] == i]
		control_samples_2 = control_samples_1[control_samples_1['chr'].isin(list(case_samples['chr'].unique()))]
		count = 0
		iters = case_samples.iterrows
		for j in range(0,args.BOOTSTRAP):
			count += 1
			if args.QUIET == False:
				if count % (args.BOOTSTRAP/10) == 0:
					print("%s %s resamples processed" % (str(j + 1),str(i)))
			stats = defaultdict(list)
			for index,row in iters():
				control_subset = filter_sample(row,control_samples_2)
				if len(control_subset) > 0:
					row = row[values]
					control_subset_2 = control_subset[values].sample()
					for k, v in dict(row-control_subset_2).items():
						stats[k].extend(v.tolist())
			vv = {"class":i,"rep":j}
			stats = dict(stats)
			for v in stats:
				vv[v] = mean(stats[v])
			df  = pd.DataFrame([vv], columns=vv.keys())
			output_file = pd.concat([output_file, df], axis =0,sort=False,ignore_index=True)

output_file.to_csv(args.OUTPUT,sep="\t",index =False)
print("--- %s seconds ---" % (time.time() - start_time))
