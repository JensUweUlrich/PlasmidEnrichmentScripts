#!/usr/bin/python

import gzip
from Bio import SeqIO
import csv
import glob
import os
import getopt,sys,re

read_ids = []


def usage():
    print("Usage:   fq_from_tsv.py -s <sequencing_summary_file> -d <fastq_dir> -o <output_file> -b <barcode>")
    print("Example: fq_from_tsv.py -s RBK403_pass_257-512.tsv -d 20220324_1028_MN24598_FAR92672_9a263022/fastq/pass -o 20220324_1028_MN24598_FAR92672_9a263022/RBK403_pass_257-512.fastq.gz -b 03")

try:
    options, remainder=getopt.getopt(sys.argv[1:], 's:d:b:o:h')

except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit()

for opt, arg in options:
    if opt in ('-s'):
        tsv_file=arg
    if opt in ('-d'):
        fastq_dir=arg
    if opt in ('-b'):
        barcode=str(arg)
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
    folder = fastq_dir + "/pass/barcode" + str(barcode)
    print(folder)
    for filename in glob.glob(folder + "/*.fastq.gz"):
        with gzip.open(filename, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if record.id in read_ids:
                    r=SeqIO.write(record, f_out, 'fastq')
                    if r!=1: print('Error while writing sequence:  ' + record.id)
                    count += 1

    for filename in glob.glob(fastq_dir + "/pass/unclassified/*.fastq.gz"):
        with gzip.open(filename, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if record.id in read_ids:
                    r=SeqIO.write(record, f_out, 'fastq')
                    if r!=1: print('Error while writing sequence:  ' + record.id)
                    count += 1
print(count)
