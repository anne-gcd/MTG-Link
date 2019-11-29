#!/usr/bin/env python3
import os
import sys
import re
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
            if re.match(self.name, record.id):
                if self.orient == "+":
                    return record.seq
                elif self.orient == "-":
                    return rc(record.seq)

    def chunk(self, c):
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
def extract_barcodes(bam, region, out_barcodes, barcodes_occ):
    command = ["BamExtractor", bam, region]
    bamextractorLog = "{}_bamextractor.log".format(datetime.now())

    with open("bam-extractor-stdout.txt", "w+") as f, open(bamextractorLog, "a") as log:
        subprocess.run(command, stdout=f, stderr=log)
        f.seek(0)

        for line in f.readlines():
            #remove the '-1' at the end of the sequence
            barcode_seq = line.split('-')[0]
            #count occurences of each barcode and save them in a dict
            if barcode_seq in barcodes_occ:
                barcodes_occ[barcode_seq] += 1
            else:
                barcodes_occ[barcode_seq] = 1

            #write the barcodes' sequences to the output file
            out_barcodes.write(barcode_seq + "\n")

    #remove the raw files obtained from BamExtractor
    subprocess.run("rm bam-extractor-stdout.txt", shell=True)
    if os.path.getsize(bamextractorLog) <= 0:
        subprocess.run(["rm", bamextractorLog])

    return out_barcodes, barcodes_occ


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
def stats_align(qry_file, ref_file, ext, prefix, out_dir):
    scriptPath = sys.path[0]
    stats_align_command = os.path.join(scriptPath, "stats_alignment.py")
    command = [stats_align_command, "-qry", qry_file, "-ref", ref_file, "-ext", ext, "-p", prefix, "-out", out_dir]
    statsLog = "{}_stats_align.log".format(datetime.now())

    with open(statsLog, "a") as log:
        subprocess.run(command, stderr=log)

    #remove the raw file obtained from statistics
    if os.path.getsize(statsLog) <= 0:
        subprocess.run(["rm", statsLog])


#----------------------------------------------------
# get_position_for_edges function
#----------------------------------------------------
def get_position_for_edges(orient1, orient2, length1, length2, k):
    #Same orientation
    if orient1 == orient2:
        beg1 = str(length1 - 2*k)
        end1 = str(length1) + "$"  
        beg2 = str(0)
        end2 = str(2*k)

    #Opposite orientation
    elif orient1 != orient2:
        beg1 = str(length1 - 2*k)
        end1 = str(length1) + "$"
        beg2 = str(length2 - 2*k)
        end2 = str(length2) + "$"

    positions = [beg1, end1, beg2, end2]
    return positions


#----------------------------------------------------
# output_gfa_with_solution function
#----------------------------------------------------
#Function to update the GFA when a solution is found for a gap
def output_gfa_with_solution(outDir, record, k, s1, s2, left_scaffold, right_scaffold, gfa_name, gfa_output_file):
    seq = record.seq
    length_seq = len(seq)
    sol_name = ".k" + str(k)
    orient = "+"

    if "bkpt2" in str(record.id):
        orient = "-"
    if "solution" in record.description:
        sol_name = record.description.split(" ")[-1] + sol_name


    sol_name = str(s1) +":"+ str(s2) + "_gf" + sol_name + orient

    #Save the found seq to a file containing all gapfill seq
    gapfill_file = gfa_name + ".gapfill_seq.fasta"
    print("Corresponding file containing all gapfill sequences: " + gapfill_file)
    with open(gapfill_file, "a") as seq_fasta:
        seq_fasta.write(">{} _ len {}".format(sol_name, length_seq))
        seq_fasta.write("\n" + str(seq) + "\n")

    with open(gfa_output_file, "a") as f:
        if length_seq < 2*k:
            print("Query length is too short (<2*k): overlap of source and destination read")
            print("Rewriting the gap line to the GFA output file...")

            #Rewrite the current G line into GFA output
            with open("tmp.gap", "r") as tmp_gap, open(gfa_output_file, "a") as f:
                out_gfa = gfapy.Gfa.from_file(gfa_output_file)
                for line in tmp_gap.readlines():
                    out_gfa.add_line(line)
                out_gfa.to_file(gfa_output_file)

        else:
            #Add the found seq (query seq) to GFA output (S line)
            out_gfa = gfapy.Gfa.from_file(gfa_output_file)
            out_gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(sol_name, length_seq, os.path.join(outDir, gapfill_file)))

            #Write the two corresponding E lines into GFA output
            pos_1 = get_position_for_edges(left_scaffold.orient, orient, left_scaffold.len, length_seq, k)
            out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(s1, sol_name, pos_1[0], pos_1[1], pos_1[2], pos_1[3]))
            pos_2 = get_position_for_edges(orient, right_scaffold.orient, length_seq, right_scaffold.len, k)
            out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(sol_name, s2, pos_2[0], pos_2[1], pos_2[2], pos_2[3]))

            out_gfa.to_file(gfa_output_file)