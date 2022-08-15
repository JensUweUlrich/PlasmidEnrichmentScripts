#!/usr/bin/python

import glob
import pysam
import csv

from ont_fast5_api.fast5_interface import get_fast5_file
from ont_fast5_api.conversion_tools.conversion_utils import get_fast5_file_list, get_progress_bar
from ont_fast5_api.fast5_file import EmptyFast5, Fast5FileTypeError
from ont_fast5_api.fast5_interface import check_file_type, MULTI_READ
from ont_fast5_api.multi_fast5 import MultiFast5File


barcodes = {}
plasmids = []
chromosomes = []

def parse_sample_sheet(data_dir):
    barcodes.clear()
    plasmids.clear()
    with open(data_dir + "/sample_sheet.csv", newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            barcodes[row['Barcode']] = row['Reference']
            for pl in row['Plasmids'].split(";"):
                plasmids.append(pl)
            for chr in row['Chromosome'].split(";"):
                chromosomes.append(chr)

def print_all_raw_data():
    fast5_filepath = "test/data/single_reads/read0.fast5" # This can be a single- or multi-read file
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            print(read.read_id, raw_data)

def get_mapping_data(data_dir):
    bam_file = data_dir + "/alignment/all.bam"
    samfile = pysam.AlignmentFile(bam_file, "rb")
    plasmid_reads = 0
    chr_reads = 0
    reads = 0
    #mapped = {}
    #mapped['plasmid'] = []
    #mapped['chromosome'] = []
    with open(data_dir + "/Plasmid_read_list.txt", 'w') as plasmid_txt:
        with open(data_dir + "/Chromosome_read_list.txt", 'w') as chr_txt:
            for read in samfile.fetch(until_eof=True):
                reads +=1
        
                if read.reference_name in plasmids:
                    plasmid_reads += 1
                    #mapped['plasmid'].append(read.query_name)
                    plasmid_txt.write(read.query_name + "\n")
                elif read.reference_name in chromosomes:
                    chr_reads += 1
                    #mapped['chromosome'].append(read.query_name)
                    chr_txt.write(read.query_name + "\n")
    print("Mapped Chromosome: " + str(chr_reads) + " / " + str(reads))
    print("Mapped Plasmid   : " + str(plasmid_reads) + " / " + str(reads))
    #return(mapped)

def main():
    data_dir = "/media/jens/INTENSO/Plasmide/Data"
    for data_folder in glob.glob(data_dir + "/2022042*"):
        parse_sample_sheet(data_folder)
        get_mapping_data(data_folder)
        #plasmid_reads = 0
        #chr_reads = 0
        #reads = 0
        #with MultiFast5File(data_folder + "/Plasmid.fast5", 'w') as plasmid_f5:
        #    with MultiFast5File(data_folder + "/Chromosome.fast5", 'w') as chr_f5:
        #        for f5_file in get_fast5_file_list(data_folder + "/fast5/", True):
                    #print(f5_file)
        ##            with get_fast5_file(f5_file, mode="r") as f5:
        #                for read in f5.get_reads():
        #                    reads +=1
        #                    #print(read.read_id)
        #                    if read.read_id in mapped['plasmid']:
        #                        #print(read.read_id)
        #                        plasmid_reads += 1
        #                        plasmid_f5.add_existing_read(read)
        #                    if read.read_id in mapped['chromosome']:
        #                        chr_reads += 1
        #                        chr_f5.add_existing_read(read)
                    #break
        #print("Mapped Chromosome: " + str(chr_reads) + " / " + str(reads))
        #print("Mapped Plasmid   : " + str(plasmid_reads) + " / " + str(reads))    
        

if __name__ == "__main__":
    main()