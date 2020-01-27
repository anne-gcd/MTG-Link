#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import re
import argparse
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="gfa2fasta.py", usage="%(prog)s -in <gfa_file> -out <output_directory>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Transform a GFA file (GFA 1.0) to a FASTA file (gaps are returned as 'Ns' regions).
                                '''))

parser.add_argument("-in", "--input", action="store", help="GFA 1.0 file (format: 'xxx.gfa')", required=True)
parser.add_argument("-out", "--outDir", action="store", help="output directory for saving the FASTA file", required=True)

args = parser.parse_args()

if re.match('^.*.gfa$', args.input) is None:
    parser.error("The suffix of the input GFA file (GFA 1.0) should be: '.gfa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
gfa_file = os.path.abspath(args.input)
if not os.path.exists(args.input):
    parser.error("The path of the input GFA file doesn't exist")
gfa_name = gfa_file.split('/')[-1]
print("\nInput GFA file: " + gfa_file)

#----------------------------------------------------
# Directory for saving results
#----------------------------------------------------
cwd = os.getcwd()
if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)
try:
    os.chdir(args.outDir)
except:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()
print("The results are saved in " + outDir)

#----------------------------------------------------
# GFA to FASTA
#----------------------------------------------------
try:
    all_seq = {}
    fasta_links = []
    assembly = ""
    i = 0
    with open(gfa_file, "r") as f:
        gfa = gfapy.Gfa.from_file(gfa_file)

        #Iterate over the 'Segment' lines of the input GFA
        for line in gfa.segments:
            seq_link = str(line).split('\t')[4].split(":")[2]

            #Open the specified seq_link only once (even if it appears for several segments in GFA)
            if seq_link not in fasta_links:
                fasta_links.append(seq_link)

        #Dictionary containing the scaffolds sequences and the gapfilled sequences
        for seq_link in fasta_links:
            with open(seq_link, "r") as input_fasta:
                current_name = ""
                current_seq = ""

                for line in input_fasta:    
                    #header line
                    if line[0] == ">":
                        current_name = (line.lstrip(">")).rstrip("\n").split(" ")[0]
                    else:
                        current_seq += line.rstrip('\n')
                        all_seq[current_name] = current_seq
                  
        #Get the paths from the 'Path' line of the input GFA
        #only the fwd gapfilled seq were kept to construct the path
        for line in gfa.paths:
            scaffolds = str(line[2]).split(',')
            overlap_lengths = str(line[3]).split(',')

            strand = scaffolds[i][-1:]
            #fwd strand
            if strand == "+":
                sequence = all_seq[scaffolds[i][:-1]]
            #rev strand
            if strand == "-":
                sequence = rc(all_seq[scaffolds[i][:-1]])

            #Attn: remove overlap only on scaffolds, not on gapfilled seq
            #initiation
            if i == 0:
                over = overlap_lengths[i][:-1]
                assembly += sequence[:-over]
            
            #last scaffold
            elif i == len(scaffolds)-1:
                over = overlap_lengths[i-1][:-1]
                assembly += sequence[over:]


            #if gapfilled sequence
            elif i % 2 == 1:
                assembly += all_seq[scaffolds[i][:-1]]

            #if scaffold
            else:
                over = overlap_lengths[i][:-1]
                assembly += (sequence[over:])[:-over]

            i += 1

    #Write the assembly sequence to the FASTA file
    fasta_name = gfa_name.split('.gfa')[0] + "_assembly.fasta"
    name = ""
    for i in range(0, len(scaffolds), 2):
        name += scaffolds[i][:-1]
        if scaffolds[i][-1:] == "+":
            name += "f_"
        else:
            name += "r_"
    name = name[:-1]    #to remove the last '_' char
    name += " len " + str(len(assembly))

    with open(fasta_name, "w") as fasta:
        fasta.write(">" + name)
        fasta.write("\n" + str(assembly) + "\n")
    

except Exception as e:
    print("\nException-")
    exc_type, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe FASTA output file is saved in " + outDir)