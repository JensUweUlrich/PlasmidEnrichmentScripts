#!/usr/bin/python

from genericpath import exists
import argparse
import shutil
import subprocess
import glob
import os
import csv
import pysam
import wget
import threading
import tarfile
import gzip
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord 


barcodes = {}
plasmids = {}
chromosomes = {}
sra = {}
min_read_length = 1000

class ThreadPool(object):
    def __init__(self):
        super(ThreadPool, self).__init__()
        self.active = []
        self.lock = threading.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
            #logging.debug('Running: %s', self.active)
    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)
            #logging.debug('Running: %s', self.active)

def execute_subcommand_thread(s, pool, cmd):
    #logging.debug('Waiting to join the pool')
    with s:
        name = threading.currentThread().getName()
        pool.makeActive(name)
        subprocess.run(cmd)
        pool.makeInactive(name)

def parse_sample_sheet(data_dir):
    barcodes.clear()
    plasmids.clear()
    with open(data_dir + "/sample_sheet.csv", newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            sra[row['Barcode']] = row['SRA']
            if row['Barcode'] == "Unclassified":
                continue
            barcodes[row['Barcode']] = row['Reference']
            for pl in row['Plasmids'].split(";"):
                plasmids[pl] = row['Reference']
            for chr in row['Chromosome'].split(";"):
                chromosomes[chr] = row['Reference']
            

def create_tsv_rd_length_plot(data_dir):
    cmd = ["Rscript", "scripts/Read_length_hist.R", "--data_dir", data_dir,"--min_read_length", str(min_read_length)]
    #print(cmd)
    subprocess.run(cmd)


def create_fastq(data_dir):
    csv_dir = data_dir + "/csv"
    fastq_dir = data_dir + "/fastq"
    #if not os.path.exists(fastq_dir):
    #    os.mkdir(fastq_dir)
    commands = []
    for filename in glob.glob(csv_dir + "/*.tsv"):
        head, tail = os.path.split(filename)
        newname = os.path.splitext(tail)[0]
        barcode = str(tail[0:6])
        cmd = ["scripts/fq_from_tsv.py", "-i", fastq_dir + "/" + barcode + ".fastq.gz", "-o", fastq_dir + "/" + newname + ".fastq.gz", "-s", filename, "-u", fastq_dir + "/Unclassified.fastq.gz"]
        #commands.append(cmd)
        #print(cmd)
        #subprocess.run(cmd)

    

    fastq_dir = fastq_dir + "/timesplit"
    if not os.path.exists(fastq_dir):
        os.mkdir(fastq_dir)

    for dir in glob.glob(csv_dir + "/timesplit/*"):
        head, tail = os.path.split(dir)
        out_dir = fastq_dir + "/" + tail
        #if int(tail) <= 960:
        #    continue
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        for filename in glob.glob(dir + "/*.tsv"):
            head, tail = os.path.split(filename)
            newname = os.path.splitext(tail)[0]
            barcode = str(tail[0:6])
            cmd = ["scripts/fq_from_tsv.py", "-i", fastq_dir + "/" + barcode + ".fastq.gz", "-o", out_dir + "/" + newname + ".fastq.gz", "-s", filename, "-u", fastq_dir + "/Unclassified.fastq.gz"]
            #print(cmd)
            #subprocess.run(cmd)
            commands.append(cmd)

    pool = ThreadPool()
    s = threading.Semaphore(8)
    i = 0
    threads = []
    for cmd in commands:
        i += 1
        t = threading.Thread(target=execute_subcommand_thread, name='thread_'+str(i), args=(s, pool, cmd))
        threads.append(t)
        t.start()
    
    for t in threads:
        t.join()

def create_enrichment_plot(data_dir):
    result_dir = data_dir + "/results/"
    cov_report_file = result_dir +  "coverage_report.csv"
    seq_summary_file = (glob.glob(data_dir + "/sequencing_summary_*.txt"))[0]
    cmd = ["Rscript", "scripts/enrichment_plot.R", "--data_dir", data_dir, "--output_dir", result_dir]
    #print(cmd)
    subprocess.run(cmd)


def create_bam(genome_dir, refname, filename, aln_dir, channels):
    mm2_cmd = ["minimap2","-ax", "map-ont", genome_dir + "/" + refname + ".fasta", filename]
    minimapoutput = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
    #samcmd = ["samtools","view", "-bh","--threads","8"]
    #samoutput = subprocess.Popen(samcmd, stdin=minimapoutput.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    bamfile = aln_dir + "/" + refname + "_" + channels + ".bam"
    samsortcmd = ["samtools", "sort","-O","BAM", "-o", bamfile]
    samsortout = subprocess.Popen(samsortcmd ,stdin=minimapoutput.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    minimapoutput.stdout.close()
    #samoutput.stdout.close()
    samsortout.communicate()
    statcommand = ["samtools", "index", bamfile]
    subprocess.run(statcommand)

def map_reads(data_dir, genome_dir):

    catcommand = ["cat", data_dir + "/fastq/*.fastq.gz", ">", data_dir + "/fastq/all.fastq.gz"]
    subprocess.run(catcommand)

    aln_dir = data_dir + "/alignment"
    fastq_dir = data_dir + "/fastq"
    result_dir = data_dir + "/results"
    if not os.path.exists(aln_dir):
        os.mkdir(aln_dir)
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)

    for filename in glob.glob(fastq_dir + "/*.fastq.gz"):
        print(filename)
        head, tail = os.path.split(filename)
        newname = os.path.splitext(tail)[0]
        newname = os.path.splitext(newname)[0]
        namesplit = newname.split("_")
        create_bam(genome_dir, barcodes[namesplit[0]], filename, aln_dir, namesplit[2])
                   
        
def create_per_read_mappings(data_dir):
    mapped_csv = {}
    sum_files = glob.glob(data_dir + "/sequencing_summary*.txt")
    with open(sum_files[0], newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
#            print(row)
            mapped_csv[row['read_id']] = {}
            mapped_csv[row['read_id']]['channel'] = row['channel']

    dir = data_dir + "/alignment"
    for filename in glob.glob(dir + "/all.bam"):
        #head, tail = os.path.split(filename)
        #newname = os.path.splitext(tail)[0]
        #ref = newname.split("_")[0]
        samfile = pysam.AlignmentFile(filename, "rb")
        plasmid_reads = 0
        chr_reads = 0
        both = 0
        none = 0
        reads = 0
        for read in samfile.fetch(until_eof=True):   
            if not 'contig' in mapped_csv[read.query_name]:
                reads +=1
                if read.reference_name in plasmids:
                    plasmid_reads += 1
                    mapped_csv[read.query_name]['contig'] = 'plasmid'
                    mapped_csv[read.query_name]['reference'] = plasmids[read.reference_name]
                elif read.reference_name in chromosomes:
                    chr_reads += 1
                    mapped_csv[read.query_name]['contig'] = 'chromosome'
                    mapped_csv[read.query_name]['reference'] = chromosomes[read.reference_name]
                else:
                    none += 1
                    mapped_csv[read.query_name]['contig'] = 'none'
                    mapped_csv[read.query_name]['reference'] = "-"
            else:
                if read.reference_name in plasmids:
                    if mapped_csv[read.query_name]['contig'] == 'chromosome':
                        plasmid_reads -= 1
                        both += 1
                        mapped_csv[read.query_name]['contig'] = 'both'
                if read.reference_name in chromosomes:
                    if mapped_csv[read.query_name]['contig'] == 'plasmid':
                        chr_reads -= 1
                        both += 1
                        mapped_csv[read.query_name]['contig'] = 'both'
                    
        print("Reads: " + str(reads))
        print("Plasmid reads: " + str(plasmid_reads))
        print("Chromosom reads: " + str(chr_reads))
        print("Both mapped: " + str(both))
        print("None mapped: " + str(none))

    with open(data_dir + "/mapping_summary.txt", 'w') as csvfile:
        csvfile.write("read_id\tchannel\treference\tcontig\n")
        #csvfile.write("read_id\tchannel\treference\tchromosome\n")
        for read_id in mapped_csv.keys():
            if 'reference' in mapped_csv[read_id]:
                csvfile.write(read_id + "\t" + str(mapped_csv[read_id]['channel']) + "\t" + mapped_csv[read_id]['reference'] + "\t" + mapped_csv[read_id]['contig'] + "\n")


def create_depth_report(data_dir, genome_dir):
    fastq_dir = data_dir + "/fastq"
    with open(data_dir + "/plasmid_depth_summary.txt", 'w') as csvfile:
        csvfile.write("reference\tchannel\ttimesplit\tmeandepth\n")
        for dir in glob.glob(fastq_dir + "/timesplit/*"):
            head, tail = os.path.split(dir)
            timepoint = tail
            print(timepoint)
            for filename in glob.glob(dir + "/*.fastq.gz"):
                head, tail = os.path.split(filename)
                newname = os.path.splitext(tail)[0]
                newname = os.path.splitext(tail)[0]
                newname = os.path.splitext(newname)[0]
                namesplit = newname.split("_")
                mm2_cmd = ["minimap2","-t","8","-ax", "map-ont", genome_dir + "/" + barcodes[namesplit[0]] + ".fasta", filename]
                minimapoutput = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
                bamfile = fastq_dir + "/" + barcodes[namesplit[0]] + "_" + namesplit[2] + ".bam"
                samsortcmd = ["samtools", "sort","-O","BAM", "-o", bamfile]
                samsortout = subprocess.Popen(samsortcmd ,stdin=minimapoutput.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                minimapoutput.stdout.close()
                samsortout.communicate()
                samcovcmd = ["samtools", "coverage", bamfile]
                samcovout = subprocess.Popen(samcovcmd ,stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
                samout,err = samcovout.communicate()
                statlines = samout.split("\n")
                meandepth = 0
                plcount = 0
                for line in statlines:
                    if str(line).startswith("#rname"):
                        continue
                    entries = str(line).split()
                    if len(entries) < 6:
                        continue
                    if entries[0] in plasmids:
                        meandepth += float(entries[6])
                        plcount += 1
                    #print(entries[0])
                    #sl = str(line).split()
                    #if len(sl) < 3:
                    #    continue
                    #if sl[3] == "mapped":
                    #    mapped = sl[0]
                os.remove(bamfile)
                meandepth /= plcount
                #print(meandepth)
                csvfile.write(barcodes[namesplit[0]] + "\t" + namesplit[2] + "\t" + str(timepoint) + "\t" + str(meandepth) + "\n")
                #return

def remove_intermediate_files(assembly_dir):
    print("Removing intermediate flye/medaka assembly files")
    for dir in glob.glob(assembly_dir + "/*0-*"):
        #print(dir)
        shutil.rmtree(dir)
    for file in glob.glob(assembly_dir + "/medaka_consensus/*.bam*"):
        #print(file)
        os.remove(file)
    for file in glob.glob(assembly_dir + "/medaka_consensus/*.hdf"):
        #print(file)
        os.remove(file)
    for file in glob.glob(assembly_dir + "/*.mmi"):
        #print(file)
        os.remove(file)
    

def assemble_reads(data_dir, genome_dir):
    assembly_dir = data_dir + "/assembly"
    if not os.path.exists(assembly_dir):
        os.mkdir(assembly_dir)

    assembly_dir = assembly_dir + "/timesplit"
    if not os.path.exists(assembly_dir):
        os.mkdir(assembly_dir)

    for dir in glob.glob(data_dir + "/fastq/timesplit/*"):
        head, tail = os.path.split(dir)
        out_dir = assembly_dir + "/" + tail
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        for filename in glob.glob(dir + "/*.fastq.gz"):
            head, tail = os.path.split(filename)
            newname = os.path.splitext(tail)[0]
            newname = os.path.splitext(newname)[0]
            namesplit = newname.split("_")
            genome_size = "100k"
            if str(namesplit[2]) == "257-512":
                genome_size = "5m"
            flye_cmd =["flye","--nano-hq",filename,"--out-dir",out_dir + "/" + newname ,"--threads","8","-g",genome_size,"--asm-coverage", "50"]
            if not os.path.exists(out_dir + "/" + newname):
                subprocess.run(flye_cmd)
                #print(flye_cmd)
                medaka_cmd = ["medaka_consensus","-i", filename,"-d", out_dir + "/" + newname + "/assembly.fasta","-o",out_dir + "/" + newname + "/medaka_consensus","-t","8","-m","r941_min_sup_g507"]
                subprocess.run(medaka_cmd)
            remove_intermediate_files(out_dir + "/" + newname)
            quast_cmd = ["quast.py","-o",out_dir + "/" + newname + "/quast_out","-r",genome_dir + "/" + barcodes[namesplit[0]] + "_plasmid.fasta","-t","8","--nanopore",filename, out_dir + "/" + newname + "/medaka_consensus/consensus.fasta"]
            if not os.path.exists(out_dir + "/" + newname + "/quast_out"):
                subprocess.run(quast_cmd)
            #mash_cmd = ["mash","dist","-i",out_dir + "/" + newname + "/medaka_consensus/consensus.fasta","/home/jens/Databases/PLSDB/plsdb.fasta.msh", ">", out_dir + "/" + newname + "/medaka_consensus/plsdb_distances.tab"]
            #pf_out_dir = out_dir + "/" + newname + "/pl_finder_out"
            #if not os.path.exists(pf_out_dir):
            #    os.mkdir(pf_out_dir)
            #plasmidfinder_cmd = ["python","/home/jens/software/plasmidfinder/plasmidfinder.py","-i",out_dir + "/" + newname + "/medaka_consensus/consensus.fasta","-o",pf_out_dir,"-p","/home/jens/software/plasmidfinder/plasmidfinder_db","-l","0.6","-t","0.95"]
            #return

def combine_quast_output(data_dir):

    with open(data_dir + "/assembly_quality_stats.txt", 'w') as outputfile:
        outputfile.write("reference\tchannel\ttimesplit\treads_mapped\tLGA50\tNGA50\tavg_coverage_depth\tref_cov1\tref_cov5\tref_cov10\tmismatches_100kb\tindels_100kb\n")
        for time_dir in glob.glob(data_dir + "/assembly/timesplit/*"):
            head, tail = os.path.split(time_dir)
            timepoint = tail
            for dir in glob.glob(time_dir + "/*"):
                head, tail = os.path.split(dir)
                namesplit = tail.split("_")
                reference = barcodes[namesplit[0]]
                channel = str(namesplit[2])
                quast_dir = dir + "/quast_out"
                stats = {}
                stats['# mismatches per 100 kbp'] = "-"
                stats['# indels per 100 kbp'] = "-"
                stats['LGA50'] = "-"
                if not os.path.exists(quast_dir + "/report.tsv"):
                    continue
                with open(quast_dir + "/report.tsv") as csvfile:
                    for line in csvfile:
                        if str(line).startswith("Assembly"):
                            continue
                        entries = str(line).rstrip("\n").split("\t")
                        stats[entries[0]] = entries[1]

                with open(quast_dir + "/reads_stats/reads_report.tsv") as csvfile:
                    for line in csvfile:
                        if str(line).startswith("Assembly"):
                            continue
                        entries = str(line).rstrip("\n").split("\t")
                        stats[entries[0]] = entries[1]
                
                #print(stats)
                outputfile.write(str(reference) + "\t" + str(channel) + "\t" + str(timepoint) + "\t" + str(stats['# reference mapped']) + "\t" + str(stats['LGA50']) + "\t" + str(stats['NGA50']) + "\t" + str(stats['Reference avg. coverage depth']) \
                    + "\t" + str(stats['Reference coverage >= 1x (%)']) + "\t" + str(stats['Reference coverage >= 5x (%)']) \
                    + "\t" + str(stats['Reference coverage >= 10x (%)']) + "\t" + str(stats['# mismatches per 100 kbp']) \
                    + "\t" + str(stats['# indels per 100 kbp']) + "\n")
                #return

def create_unblocked_plasmid_fastq_files(data_dir):
    unblocked_read_ids = []
    lengths = 0
    unblock_reads = 0
    with open(data_dir + "/other_reports/read_until_decision_stats.csv", newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=';')
        for row in reader:
            if row['decision'] == "unblock":
                if not row['read_id'] in unblocked_read_ids:
                    lengths += int(row["sequence_length"])
                    unblock_reads += 1
                    unblocked_read_ids.append(row['read_id'])
    avg_len = float(lengths) / float(unblock_reads)
    print("avg. unblock length: " + str(avg_len))

    dir = data_dir + "/alignment"
    unblocked_plasmid_read_ids = {}
    for filename in glob.glob(dir + "/all.bam"):
        samfile = pysam.AlignmentFile(filename, "rb")
        for read in samfile.fetch(until_eof=True):
            if read.reference_name in plasmids:
                refname = plasmids[read.reference_name]
                if not refname in unblocked_plasmid_read_ids:
                    unblocked_plasmid_read_ids[refname] = []
                if read.query_name in unblocked_read_ids:
                    unblocked_plasmid_read_ids[refname].append(read.query_name)

    for refname in unblocked_plasmid_read_ids:
        with gzip.open(data_dir + "/fastq/" + refname + "_unblocked_plasmid.fastq.gz", 'wt') as f_out:            
            with gzip.open(data_dir + "/fastq/all.fastq.gz", "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    if record.id in unblocked_plasmid_read_ids[refname]:
                        r=SeqIO.write(record, f_out, 'fastq')
                        if r!=1: print('Error while writing sequence:  ' + record.id)


#def combine_experiment_plots(plasmid_dir):
#    result_dir = plasmid_dir + "/results/"
#    cov_report_file = result_dir +  "coverage_report.csv"
#    seq_summary_file = (glob.glob(data_dir + "/sequencing_summary_*.txt"))[0]
#    cmd = ["Rscript", "scripts/enrichment_plot.R", "--data_dir", data_dir, "--output_dir", result_dir]
#    print(cmd)
#    subprocess.run(cmd)

def map_unblocked_plasmid2chromosome(data_dir, genome_dir):

    with open(data_dir + "/unblocked_plasmid_chromosome_regions.txt", 'w') as outputfile:
        outputfile.write("reference\tcontig\tstart\tend\n")
        for bc in barcodes.keys():
            mm2_cmd = ["minimap2","-ax", "map-ont", genome_dir + "/" + barcodes[bc] + "_chromosome.fasta", data_dir + "/fastq/" + barcodes[bc] + "_unblocked_plasmid.fastq.gz"]
            minimapoutput = subprocess.Popen(mm2_cmd, stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
            bamfile = data_dir + "/alignment/unblocked_plasmid_" + barcodes[bc] + ".bam"
            samsortcmd = ["samtools", "sort","-O","BAM", "-o", bamfile]
            samsortout = subprocess.Popen(samsortcmd ,stdin=minimapoutput.stdout, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            minimapoutput.stdout.close()
            samsortout.communicate()
            statcommand = ["samtools", "index", bamfile]
            subprocess.run(statcommand)
            samdepthcmd = ["samtools", "depth", bamfile]
            samdepthout = subprocess.Popen(samdepthcmd ,stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, universal_newlines=True)
            samout,err = samdepthout.communicate()
            statlines = samout.split("\n")
            start = 0
            end = 0
            regions = {}
            count10 = 0
            ref = ""
            for line in statlines:
                #print(line)
                entries = str(line).split()
                if len(entries) < 3:
                    continue
                if int(entries[2]) >= 10:
                    count10 = 0
                    if start == 0:
                        start = int(entries[1])
                        ref = entries[0]
                        continue
                    else:
                        if ref != entries[0]:
                            if not ref in regions:
                                regions[ref] = []
                            regions[ref].append(str(start) + "-" + str(end))
                            count10 = 0
                            start = int(entries[1])
                            ref = entries[0]
                            end = 0
                        else:
                            end = int(entries[1])
                        continue
                else:
                    if start > 0 and end > 0:
                        count10 += 1
                        if count10 >= 10:
                            if not ref in regions:
                                regions[ref] = []
                            regions[ref].append(str(start) + "-" + str(end))
                            count10 = 0
                            start = 0
                            end = 0
            print(barcodes[bc])
            l = 0
            
            for ref in regions:
                for reg in regions[ref]:
                    start, end = str(reg).split("-")
                    outputfile.write(barcodes[bc] + "\t" +  ref + "\t" + start + "\t" + end + "\n")

            regionnr = 0
            with open(data_dir + "/fastq/" + barcodes[bc] +"_unblocked_plasmid_chr_regions.fasta", 'wt') as f_out:
                with open(genome_dir + "/" + barcodes[bc] + "_chromosome.fasta", "rt") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        for ref in regions:
                            if record.id == ref:
                                for reg in regions[ref]:
                                    start, end = str(reg).split("-")
                                    #print(record.seq[int(start)-1 : int(end)])
                                    regionnr += 1
                                    new_record = SeqIO.SeqRecord(Seq(record.seq[int(start)-1 : int(end)]), id="region_"+str(regionnr))
                                    r=SeqIO.write(new_record, f_out, 'fasta')
                                    if r!=1: print('Error while writing sequence:  ' + new_record.id)

def prepare_data(dir):

    for zipped in glob.glob(dir + "/*.tar.gz"):                
        file = tarfile.open(zipped)
        file.extractall(dir + "/.")
        file.close()

    parse_sample_sheet(dir)

    fastq_dir = dir + "/fastq"
    if not os.path.exists(fastq_dir):
        os.mkdir(fastq_dir)

    for bc in sra.keys():
        if bc == "RBK401" or bc == "RBK402":
            continue
        cmd = ["prefetch",sra[bc]]
        subprocess.run(cmd)
        cmd = ["fastq-dump", "--outdir", fastq_dir, "--gzip", sra[bc] + "/" + sra[bc] + ".sra"]
        subprocess.run(cmd)
        os.rename(fastq_dir + "/" + sra[bc] + ".fastq.gz", fastq_dir + "/" + bc + ".fastq.gz")


def storeRefSeqFasta(output_file, records):
    with open(output_file, 'wt') as f_out:
        for rec in records:
            r=SeqIO.write(rec, f_out, 'fasta')
            if r!=1: print('Error while writing sequence:  ' + rec.id)

def get_reference_genomes(dir):

    orgs = {}

    with open(dir + "/references.csv", newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        for row in reader:
            r = {}
            r['name'] = row['name']
            r['chromosome'] = row['chromosome'].split(";")
            r['plasmids'] = row['plasmids'].split(";")
            orgs[row['ID']] = r


    genome_dir = dir + "/Genomes"
    if not os.path.exists(genome_dir):
        os.mkdir(genome_dir)

    Entrez.email = "jensenuk83@gmail.com"
    for sp in orgs:
        # create chromosome ref seq file
        records = []
        for ch in orgs[sp]['chromosome']:
            handle = Entrez.efetch(db="nucleotide", id=ch,rettype="fasta", retmode="text")
            for record in list(SeqIO.parse(handle, "fasta")):
                records.append(record)
            handle.close()

        storeRefSeqFasta(genome_dir + "/" + orgs[sp]['name'] + "_chromosome.fasta", records)
        plasmids = []
        # create plasmid ref seq file
        for pl in orgs[sp]['plasmids']:
            handle = Entrez.efetch(db="nucleotide", id=pl,rettype="fasta", retmode="text")
            plasmids.append(SeqIO.read(handle, "fasta"))
            handle.close()
        storeRefSeqFasta(genome_dir + "/" + orgs[sp]['name'] + "_plasmid.fasta", plasmids)

        for r in plasmids:
            records.append(r)

        storeRefSeqFasta(genome_dir + "/" + orgs[sp]['name'] + ".fasta", records)
    return genome_dir

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", help="Data directory", type=str)
    args = parser.parse_args()
 
    genome_dir = get_reference_genomes(args.data_dir)
    data_dir = args.data_dir
    for data_folder in glob.glob(data_dir + "/20220321*"):
        print(data_folder)
        prepare_data(data_folder)
        create_tsv_rd_length_plot(data_folder)
        create_fastq(data_folder)
        return
        map_reads(data_folder, genome_dir)
        create_per_read_mappings(data_folder)
        create_unblocked_plasmid_fastq_files(data_folder)
        map_unblocked_plasmid2chromosome(data_folder, genome_dir)
        create_depth_report(data_folder, genome_dir)
        create_enrichment_plot(data_folder)
        assemble_reads(data_folder, genome_dir)
        combine_quast_output(data_folder)
#    combine_experiment_plots("media/jens/INTENSO/Plasmide")

if __name__ == "__main__":
    main()
