#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import re
import argparse
import gfapy
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="fasta2gfa.py", usage="%(prog)s -in <fasta_file> -out <output_directory> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Transform a FASTA file with sequences containing 'Ns' regions to a GFA file ('Ns' regions are treated as gaps).
                                We can filter the 'Ns' regions by their size (e.g. gap sizes) and by the contigs' sizes on both sides (long enough for ex to get enough barcodes) by providing values for the following arguments:
                                    - '-min': minimum size gap
                                    - '-max': maximum size gao
                                    - '-contigs': minimum size of the contigs on both sides of the Ns regions
                                '''))

parser.add_argument("-in", "--input", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta')", required=True)
parser.add_argument("-min", "--min", action="store", type=int, help="minimum size of the Ns region to treat/process as a gap")
parser.add_argument("-max", "--max", action="store", type=int, help="maximum size of the Ns region to treat/process as a gap")
parser.add_argument("-contigs", "--contigs_size", action="store", type=int, help="minimum size of the contigs on both sides of the Ns region to treat/process the Ns region as a gap")
parser.add_argument("-out", "--outDir", action="store", help="output directory for saving the GFA file", required=True)

args = parser.parse_args()

if re.match('^.*.fasta$', args.input) is None:
    parser.error("The suffix of the input FASTA file should be: '.fasta'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
fasta_file = os.path.abspath(args.input)
if not os.path.exists(args.input):
    parser.error("The path of the input FASTA file doesn't exist")
fasta_name = fasta_file.split('/')[-1]
print("\nInput FASTA file: " + fasta_file)

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
# FASTA to GFA
#----------------------------------------------------
try:
    #Iterate over the scaffolds in the FASTA file
    with open(fasta_file, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            name = record.id
            seq = record.seq  

            gap_count = 0
            gaps = []
            start_index_gap = []
            end_index_gap = []

            #If Ns in the sequence, get the gaps' sequences
            end = 0
            while "N" in str(seq[end:]):
                match = re.search("N+", str(seq[end:]))
                gap = match.group()
                start = match.start() + end
                end = match.end() + end

                gaps.append(gap)
                start_index_gap.append(start)
                end_index_gap.append(end)

            #If 'min' < length of gap < 'max' AND left and right sequences > 'contigs_size', write the left and right sequences (on both sides of the gap) to the GFA file, and to FASTA files
            i = 0
            for gap in gaps:

                #get the left and right sequences (on both sides of the gap)
                if len(gaps) == 1:          #only one gap in the scaffold
                    left_seq = str(seq[:start_index_gap[0]])
                    right_seq = str(seq[end_index_gap[0]:])
                elif i == 0:                  #first gap
                    left_seq = str(seq[:start_index_gap[i]])
                    right_seq = str(seq[end_index_gap[i]:start_index_gap[i+1]])
                elif i == (len(gaps) - 1):  #last gap
                    left_seq = right_seq            #current left_seq = previous right_seq
                    right_seq = str(seq[end_index_gap[i]:])
                else:
                    left_seq = right_seq            #current left_seq = previous right_seq
                    right_seq = str(seq[end_index_gap[i]:start_index_gap[i+1]])

                i+= 1

                #check the sequences meet the required conditions
                if len(gap) >= args.min and len(gap) <= args.max and len(left_seq) > args.contigs_size and len(right_seq) > args.contigs_size:
                    gap_count += 1
                    os.chdir(outDir)

                    #----------------------------------------------------
                    # FASTA files (left/right sequences)
                    #----------------------------------------------------
                    #Directory for saving the corresponding FASTA files
                    out_fasta_file = os.path.abspath(fasta_name.split(".fasta")[0] + "_scaff_" + str(name) + "_gaps_" + str(args.min) + ":" + str(args.max) + "_contigs_" + str(args.contigs_size) + ".fasta")
                    print("\nOutput FASTA file: " + out_fasta_file)
                                
                    #Save the left and right sequences to FASTA files
                    left_name = str(name) + "_gap" + str(gap_count) + "-L"
                    left_contig = left_name + "+"
                    len_left = len(left_seq)

                    right_name = str(name) + "_gap" + str(gap_count) + "-R"
                    right_contig = right_name + "+"
                    len_right = len(right_seq)

                    with open(out_fasta_file, "a") as output:
                        output.write(">{} _ len {}".format(left_name, len_left))
                        output.write("\n" + str(left_seq) + "\n")
                        output.write(">{} _ len {}".format(right_name, len_right))
                        output.write("\n" + str(right_seq) + "\n")

                    #----------------------------------------------------
                    # GFA file
                    #----------------------------------------------------
                    os.chdir(outDir)

                    gfa_file = os.path.abspath(fasta_name.split(".fasta")[0] + "_scaff_" + str(name) + "_gaps_" + str(args.min) + ":" + str(args.max) + "_contigs_" + str(args.contigs_size) + ".gfa")
                    print("\nOutput GFA file: " + gfa_file)

                    #Initiate GFA file
                    if not os.path.exists(gfa_file):
                        with open(gfa_file, "w") as f:
                            gfa = gfapy.Gfa()
                            gfa.add_line("H\tVN:Z:2.0")
                            gfa.to_file(gfa_file)

                    #Add corresponding lines to GFA file
                    with open(gfa_file, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfa_file)
                        gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(left_name, len_left, os.path.join(outDir, out_fasta_file)))
                        gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(right_name, len_right, os.path.join(outDir, out_fasta_file)))
                        gfa.add_line("G\t*\t{}\t{}\t{}\t*".format(left_contig, right_contig, len(gap)))
                        gfa.to_file(gfa_file)


except Exception as e:
    print("\nException-")
    exc_type, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe GFA output files and the corresponding FASTA files are saved in " + outDir)