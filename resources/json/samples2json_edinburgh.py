#!/apps/languages/anaconda3/bin/python


import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--fastq_dir", help="Required. the FULL path to the fastq folder")
args = parser.parse_args()

assert args.fastq_dir is not None, "please provide the path to the fastq folder"


## default dictionary is quite useful!

FILES = defaultdict(lambda: defaultdict(list))

## build the dictionary with full path for each fastq.gz file
for root, dirs, files in os.walk(args.fastq_dir):
	for file in files:
		if file.endswith("fastq.gz"):
			full_path = join(root, file)
			#R1 will be forward reads, R2 will be reverse reads
			m = re.search(r"([0-9]{5})_.+_([0-9]{1})_.+([12]).fastq.gz", file)
# re.search(r"-([A-Z]{4})-([A-Z0-9a-z]{4})-[0-9A-Z]{2}-.+_(L[0-9]{3})_(R[12])_[0-9]{3}.fastq.gz", file)
#17115_150609_D00200_0263_BC7187ANXX_1_IL-TP-012_1.sanfastq.gz
			if m:
                                sample = m.group(1)
                                lane = m.group(2)
                                # R1 or R2 for forward and reverse read
                                reads = m.group(3)  
                                FILES[sample][reads].append(full_path)

# make sure R1 and R2 from different lanes are ordered in the same way
# e.g. L001_R1 should pair with L001_R2, L002_R1 pair with L002_R2		

FILES_sorted = defaultdict(lambda: defaultdict(list))

for sample in FILES.keys():
		for read in FILES[sample]:
			FILES_sorted[sample][read] = sorted(FILES[sample][read])



print()
print ("total {} unique samples will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
for sample in FILES_sorted.keys():
	for read in sorted(FILES_sorted[sample]):
		print ("{sample} {read} has {n} fastq".format(sample = sample, read = read, n = len(FILES_sorted[sample][read])))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES_sorted, indent = 4, sort_keys=True)
open('samples_edinburgh.json', 'w').writelines(js)
