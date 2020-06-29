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
parser = argparse.ArgumentParser(prog="paths2gfa.py", usage="%(prog)s -in <fasta_file> -paths <paths_file> -out <output_directory>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Transform a file containing the paths between scaffolds to a GFA file
                                '''))

parser.add_argument("-in", dest="input", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta')", required=True)
parser.add_argument("-paths", dest="paths", action="store", help="File containing the paths between scaffolds (obtained from the matrix) (format: 'xxx.paths.txt')", required=True)
parser.add_argument("-out", dest="outDir", action="store", help="Output directory for saving the GFA file and the corresponding FASTA file", required=True)

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

paths_file = os.path.abspath(args.paths)
paths_name = (paths_file.split("/")[-1]).split("paths")[0]
if not os.path.exists(args.paths):
        parser.error("The path of the input paths' file doesn't exist")
print("Input paths' file: " + paths_file)

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
print("\nThe results are saved in " + outDir)

#----------------------------------------------------
# PATHS to GFA
#----------------------------------------------------
try:
    out_fasta = outDir + "/" + paths_name + "scaffolds.fasta"
    gfa_file = outDir + "/" + paths_name + "gfa"

    #Iterate over the scaffolds in the PATH
    with open(paths_file, "r") as paths:
        for line in paths:
            line = line.split("*")[-1]
            list_scaff = line.split("+")

            for scaff_orient in list_scaff:
                scaffold = scaff_orient.split("(")[0]

                #Find sequence of scaffold in FASTA file
                with open(fasta_file, "r") as fasta:
                    for record in SeqIO.parse(fasta, "fasta"):
                        if re.match(str(scaffold), record.id):
                            seq = record.seq

                #----------------------------------------------------
                # FASTA file
                #----------------------------------------------------            
                print("\nOutput FASTA file: " + out_fasta)

                #Save the scaffolds' sequence to FASTA file
                with open(out_fasta, "a") as fasta_out:
                    fasta_out.write(">{} _ len {}".format(scaffold, len(seq)))
                    fasta_out.write("\n" + str(seq) + "\n")

                #----------------------------------------------------
                # GFA file - Segment lines
                #----------------------------------------------------
                print("Output GFA file: " + gfa_file)

                #Initiate GFA file
                if not os.path.exists(gfa_file):
                    with open(gfa_file, "w") as f:
                        gfa = gfapy.Gfa()
                        gfa.add_line("H\tVN:Z:2.0")
                        gfa.to_file(gfa_file)

                #Add corresponding S lines (Segments) to GFA file, if not already in GFA file
                with open(gfa_file, "a") as f:
                    gfa = gfapy.Gfa.from_file(gfa_file)
                    if str(scaffold) not in gfa.segments:
                        gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(str(scaffold), str(len(seq)), os.path.join(outDir, out_fasta)))                  
                        gfa.to_file(gfa_file)

            #----------------------------------------------------
            # GFA file - Gap lines
            #----------------------------------------------------
            #Iterate over the links between scaffolds, in the PATH
            for i in range(len(list_scaff)-1):

                #Orientations scaffolds Left and Right
                orient_L = re.search("[fr]", str(list_scaff[i]))
                orient_L = orient_L.group()
                if orient_L == "f":
                    left_orient = "+"
                elif orient_L == "r":
                    left_orient = "-"

                orient_R = re.search("[fr]", str(list_scaff[i+1]))
                orient_R = orient_R.group()
                if orient_R == "f":
       	       	    right_orient = "+"
       	       	elif orient_R == "r":
       	       	    right_orient = "-"

                #Name scaffolds Left and Right
                left_scaff = list_scaff[i].split("(")[0] + left_orient
                right_scaff = list_scaff[i+1].split("(")[0] + right_orient

                #Add corresponding G lines (Gaps) to GFA file (dist=0 means len(gap) unknown)
                with open(gfa_file, "a") as f:
                    gfa = gfapy.Gfa.from_file(gfa_file)
                    gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(left_scaff), str(right_scaff)))
                    gfa.to_file(gfa_file)


except Exception as e:
    print("\nException-")
    exc_type, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    sys.exit(1)