#!/usr/bin/python

import gzip
from Bio import SeqIO
import csv
import glob
import os
import getopt,sys,re

read_ids = []


def usage():
    print("Usage:   fq_from_tsv.py -s <sequencing_summary_file> -i <input_fastq.gz_file> -o <output_file> -u <unclassified_file>")
    print("Example: fq_from_tsv.py -s RBK403_pass_257-512.tsv -i fastq/RBK403.fastq.gz -o fastq/RBK403_pass_257-512.fastq.gz -u Unclassified.fastq.gz")

try:
    options, remainder=getopt.getopt(sys.argv[1:], 's:d:b:o:h')

except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit()

for opt, arg in options:
    if opt in ('-s'):
        tsv_file=arg
    if opt in ('-i'):
        input_file=arg
    if opt in ('-u'):
        unclassified=str(arg)
    if opt in ('-h'):
        usage()
        sys.exit()
    elif opt in ('-o'):
        output_file=arg

with open(tsv_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter='\t')
    for row in reader:
        read_ids.append(row['read_id'])

print(len(read_ids))
count = 0
with gzip.open(output_file, 'wt') as f_out:
    #folder = fastq_dir + "/pass/barcode" + str(barcode)
    #print(folder)
    #for filename in glob.glob(folder + "/*.fastq.gz"):
    with gzip.open(input_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in read_ids:
                r=SeqIO.write(record, f_out, 'fastq')
                if r!=1: print('Error while writing sequence:  ' + record.id)
                count += 1

    #for filename in glob.glob(fastq_dir + "/pass/unclassified/*.fastq.gz"):
    with gzip.open(unclassified, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in read_ids:
                r=SeqIO.write(record, f_out, 'fastq')
                if r!=1: print('Error while writing sequence:  ' + record.id)
                count += 1
print(count)
