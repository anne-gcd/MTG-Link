import os
import sys
import subprocess
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO
from datetime import datetime


#----------------------------------------------------
# Gap class
#----------------------------------------------------
class Gap:
    def __init__(self, gap):
        self.id = gap.gid
        self.length = gap.disp
        self.left = gap.sid1
        self.right = gap.sid2

    def label(self):
        if self.id == "*":
            return str(self.left) +"_"+ str(self.right)
        else:
            return str(self.id)

    def info(self):
        if self.id == "*":
            print("\nWORKING ON GAP: between scaffolds {} & {}; length {}\n".format(self.left, self.right, self.length))
        else:
            print("\nWORKING ON GAP: {}; length {}\n".format(self.id, self.length))


#----------------------------------------------------
# Scaffold class
#----------------------------------------------------
class Scaffold(Gap):
    def __init__(self, gap, scaffold):
        super().__init__(gap)
        self.gap = gap
        self.scaffold = scaffold
        self.name = scaffold.name
        self.orient = scaffold.orient
        self.len = scaffold.line.slen
        self.seq_path = scaffold.line.UR
    
    def sequence(self):
        for record in SeqIO.parse(self.seq_path, "fasta"):
            if self.orient == "+":
                return record.seq
            elif self.orient == "-":
                return rc(record.seq)

    def chunk(self, c):
        if c > self.len:
            sys.stderr.write("Error: The chunk size must be smaller than the length of the scaffolds\n")
            sys.exit(2)
        else:
            if (self.orient == "+" and self.scaffold == self.left) or (self.orient == "-" and self.scaffold == self.right):   #if left_fwd or right_rev
                start = self.len - c
                end = self.len
            elif (self.orient == "+" and self.scaffold == self.right) or (self.orient == "-" and self.scaffold == self.left):  #if right_fwd or left_rev
                start = 0
                end = c
            return str(self.name) +":"+ str(start) +"-"+ str(end)




#----------------------------------------------------
# extract_barcodes function
#----------------------------------------------------
#Function to extract the barcodes of reads mapping on chunks, with BamExtractor 
def extract_barcodes(bam, region, freq, out_barcodes):
    command = ["BamExtractor", bam, region]
    bamextractorLog = "{}_bamextractor.log".format(datetime.now())
    barcodes_occ = {}

    with open("bam-extractor-stdout.txt", "w+") as f, open(bamextractorLog, "a") as log:
        subprocess.run(command, stdout=f, stderr=log)
        f.seek(0)

        for line in f.readlines():
            #remove the '-1' at the end of the sequence
            barcode_seq = line.split('-')[0]
            #count occurences of each barcode
            if barcode_seq in barcodes_occ:
                barcodes_occ[barcode_seq] += 1
            else:
                barcodes_occ[barcode_seq] = 1
        
    #filter barcodes by freq
    for (barcode, occurences) in barcodes_occ.items():
        if occurences >= freq:
            out_barcodes.write(barcode + "\n")

    #remove the raw files obtained from BamExtractor
    subprocess.run("rm bam-extractor-stdout.txt", shell=True)
    if os.path.getsize(bamextractorLog) <= 0:
        subprocess.run(["rm", bamextractorLog])

    return out_barcodes


#----------------------------------------------------
# get_reads function
#----------------------------------------------------
#Function to extract the reads associated with the barcodes
def get_reads(reads, index, barcodes, out_reads):    
    command = ["GetReads", "-reads", reads, "-index", index, "-barcodes", barcodes]
    getreadsLog = "{}_getreads.log".format(datetime.now())

    with open(getreadsLog, "a") as log:
        subprocess.run(command, stdout=out_reads, stderr=log)

    #remove the raw file obtained from GetReads
    if os.path.getsize(getreadsLog) <= 0:
        subprocess.run(["rm", getreadsLog])

    return out_reads


#----------------------------------------------------
# mtg_fill function
#----------------------------------------------------
#Function to execute MindTheGap fill module
def mtg_fill(input_file, bkpt, k, a, max_nodes, max_length, nb_cores, max_memory, verbose, output_prefix):
    command = ["MindTheGap", "fill", "-in", input_file, "-bkpt", bkpt, "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(max_nodes), "-max-length", str(max_length), \
               "-nb-cores", str(nb_cores), "-max-memory", str(max_memory), "-verbose", str(verbose), "-out", output_prefix]
    mtgfillLog = "{}_mtgfill.log".format(datetime.now())

    with open(mtgfillLog, "a") as log:
        subprocess.run(command, stderr=log)
        output = subprocess.check_output(command)

    #remove the raw files obtained from MindTheGap
    subprocess.run("rm -f *.h5", shell=True)
    if os.path.getsize(mtgfillLog) <= 0:
        subprocess.run(["rm", mtgfillLog])

    return output


#----------------------------------------------------
# stats_align function
#----------------------------------------------------
#Function to do statistics on the alignment of a reference sequence and query sequences
def stats_align(qry_file, ref_file, prefix, out_dir):
    scriptPath = sys.path[0]
    stats_align_command = os.path.join(scriptPath, "stats_alignment.py")
    command = [stats_align_command, "-qry", qry_file, "-ref", ref_file, "-p", prefix, "-out", out_dir]
    statsLog = "{}_stats_align.log".format(datetime.now())

    with open(statsLog, "a") as log:
        subprocess.run(command, stderr=log)

    #remove the raw file obtained from statistics
    if os.path.getsize(statsLog) <= 0:
        subprocess.run(["rm", statsLog])