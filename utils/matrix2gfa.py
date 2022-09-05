#!/usr/bin/env python3
#*****************************************************************************
#  Name: MTG-Link
#  Description: Local assembly tool for linked-reads data
#  Copyright (C) 2020 INRAE
#  Author: Anne Guichard
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************

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
parser = argparse.ArgumentParser(prog="matrix2gfa.py", usage="%(prog)s -in <fastaFile.fasta> -matrix <matrixFile.matrix> -out <outDir> -threshold <THRESHOLD>",\
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a file containing the matrix (links between the ends of the scaffolds/contigs) to a GFA file"))

parser.add_argument("-fa", dest="fasta", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-matrix", dest="matrix", action="store", help="File containing the links between the ends of the scaffolds/contigs in tabular format (matrix)", required=True)
parser.add_argument("-threshold", dest="threshold", type=int,  action="store", help="Minimal number of links two scaffolds must share to try to fill the gap between them", required=False, default=10)
parser.add_argument("-out", dest="outDir", action="store", help="Directory to save the output GFA file and gap flanking sequences FASTA files", required=True)

args = parser.parse_args()

if (re.match('^.*.fasta$', args.fasta) is None) and (re.match('^.*.fa$', args.fasta) is None):
    parser.error("The suffix of the input FASTA file should be: '.fasta' or '.fa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
# FASTA file containing the sequences of the scaffolds obtained from the assembly.
fastaFile = os.path.abspath(args.fasta)
if not os.path.exists(args.fasta):
    parser.error("The path of the input FASTA file doesn't exist")
fasta_name = fastaFile.split('/')[-1]
print("\nInput FASTA file: " + fastaFile)
fastaDict= SeqIO.index(fastaFile, "fasta")

# Matrix file containing the links between the ends of the contigs in tabular format.
matrixFile = os.path.abspath(args.matrix)
matrix_name = (matrixFile.split("/")[-1]).split(".matrix")[0]
if not os.path.exists(args.matrix):
        parser.error("The path of the input matrix file doesn't exist")
print("Input matrix file: " + matrixFile)

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
# MATRIX to GFA
#----------------------------------------------------
stored_scaff={}
start_re = re.compile('0-\d+')

try:
    # Iterate over the links in the matrix file. 
    with open(matrixFile, "r") as mat:
        for line in mat:
            (scaff1,scaff2,links) = line.split(" ")

            # Keep only the links over a certain threshold.
            if (int(links) >= args.threshold):
                (scaff1_name, scaff1_pos) = scaff1.split(":")
                (scaff2_name, scaff2_pos) = scaff2.split(":")

                # If both scaffolds are the same, do not consider this link.
                if (scaff1_name == scaff2_name):
                    continue

                # Find the orientations of both scaffolds.
                ##If len(scaff) too short, the orientation of the scaffold is unknown, it can be 'both' + and -
                scaff1_orient = '+'
                scaff2_orient = '-'
                if (start_re.match(scaff1_pos)):
                    if (scaff1_pos == "0-" + str(len(fastaDict[scaff1_name]))):
                        scaff1_orient = 'both'
                    else: 
                        scaff1_orient = '-'
                if (start_re.match(scaff2_pos)):
                    if (scaff2_pos == "0-" + str(len(fastaDict[scaff2_name]))):
                         scaff2_orient = 'both'
                    else: 
                        scaff2_orient = '+'

                leftScaffold = str(scaff1_name) + str(scaff1_orient)
                rightScaffold = str(scaff2_name) + str(scaff2_orient)

                #----------------------------------------------------
                # FASTA file (left/right sequences)
                #----------------------------------------------------
                # Output FASTA files containing the flanking sequences of the current gap.
                ## Left flanking scaffold.
                if (not (scaff1_name in stored_scaff)):
                    leftSeq = fastaDict[scaff1_name]
                    left_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_" + str(scaff1_name) + ".scaffold.fasta")
                    left = open(left_fastaFile, "w")
                    leftSeq.description = "_ len " + str(len(leftSeq))
                    SeqIO.write(leftSeq, left, "fasta")
                    stored_scaff[scaff1_name]=1

                ## Right flanking scaffold.
                if (not (scaff2_name in stored_scaff)):
                    rightSeq = fastaDict[scaff2_name]
                    right_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_" + str(scaff2_name) + ".scaffold.fasta")
                    right = open(right_fastaFile, "w")
                    rightSeq.description = "_ len " + str(len(rightSeq))
                    SeqIO.write(rightSeq, right, "fasta")
                    stored_scaff[scaff2_name]=1  
                
                #----------------------------------------------------
                # GFA file
                #----------------------------------------------------
                os.chdir(outDir)
                gfaFile = os.path.abspath(fasta_name.split(".fa")[0] +"_"+ matrixFile.split("/")[-1].split(".matrix")[0] + "_threshold_" + str(args.threshold) + ".gfa")

                # Initiate the GFA file.
                if not os.path.exists(gfaFile):
                    with open(gfaFile, "w") as f:
                        gfa = gfapy.Gfa()
                        gfa.add_line("H\tVN:Z:2.0")
                        gfa.to_file(gfaFile)

                # Add corresponding lines to the GFA file.
                with open(gfaFile, "a") as f:
                    gfa = gfapy.Gfa.from_file(gfaFile)

                    # S lines (Segments) of the GFA file.
                    ##Left
                    if str(scaff1_name) not in gfa.segments:
                        gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(scaff1_name, str(len(leftSeq)), os.path.join(outDir, left_fastaFile)))
                    ##Right
                    if str(scaff2_name) not in gfa.segments:
                        gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(scaff2_name, str(len(rightSeq)), os.path.join(outDir, right_fastaFile)))

                    # G lines (Gaps) of the GFA file.
                    ##Orientation of scaff1 unknown (it can be 'both' + and -)
                    if (scaff1_orient == 'both'):
                        if (scaff2_orient == 'both'):
                            gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(scaff1_name)+"+", str(scaff2_name)+"+"))
                            gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(scaff1_name)+"+", str(scaff2_name)+"-"))
                            gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(scaff1_name)+"-", str(scaff2_name)+"+"))
                            gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(scaff1_name)+"-", str(scaff2_name)+"-"))
                        else :
                            gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(scaff1_name)+"+", str(rightScaffold)))
                            gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(scaff1_name)+"-", str(rightScaffold)))

                    ##Orientation of scaff2 unknown (it can be 'both' + and -)
                    elif (scaff2_orient == 'both'):
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(scaff2_name)+"+"))
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(scaff2_name)+"-"))

                    ##Orientations of both scaffolds are known
                    else:
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))

                    gfa.to_file(gfaFile)


except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe output GFA file and the corresponding FASTA files are saved in " + outDir)
print("Output GFA file: " + gfaFile)

