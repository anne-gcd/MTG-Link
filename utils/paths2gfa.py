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
parser = argparse.ArgumentParser(prog="paths2gfa.py", usage="%(prog)s -fa <fastaFile.fasta> -paths <pathsFile.txt> -out <outDir>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a file containing the paths between scaffolds to a GFA file"))

parser.add_argument("-fa", dest="fasta", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-paths", dest="paths", action="store", help="File containing the paths between scaffolds (obtained from the matrix) (format: 'xxx.paths.txt')", required=True)
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

# File containing the paths between scaffolds.
pathsFile = os.path.abspath(args.paths)
paths_name = (pathsFile.split("/")[-1]).split("paths")[0]
if not os.path.exists(args.paths):
        parser.error("The path of the input paths' file doesn't exist")
print("Input PATHS file: " + pathsFile)

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
stored_scaff={}

try:
    # Iterate over the lines of the PATHS file and get a list of the scaffolds.
    with open(pathsFile, "r") as paths:
        for line in paths:
            line = line.split("*")[-1]
            scaffoldsList = line.split("+")

            # Iterate over the scaffolds in the 'scaffoldsList' list. 
            for scaffold in scaffoldsList:
                scaffoldID = scaffold.split("(")[0]

                # Find the sequence of the current scaffold in the input FASTA file.
                with open(fastaFile, "r") as fasta:
                    for record in SeqIO.parse(fasta, "fasta"):
                        if re.match(str(scaffoldID), record.id):
                            seq = record.seq
                            break

                #----------------------------------------------------
                # FASTA file
                #----------------------------------------------------
                # Output FASTA file containing the sequence of the current scaffold.
                if (not (scaffoldID in stored_scaff)):
                    scaffold_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] +"_"+ str(record.id) +".scaffold.fasta")
                    with open(scaffold_fastaFile, "w") as scaffold_fasta:
                        scaffold_fasta.write(">{} _ len {}".format(scaffoldID, len(seq)))
                        scaffold_fasta.write("\n" + str(seq) + "\n")
                    stored_scaff[scaffoldID]=1

                #----------------------------------------------------
                # GFA file - Segment lines
                #----------------------------------------------------
                os.chdir(outDir)
                gfaFile = os.path.abspath(fasta_name.split(".fa")[0] +"_"+ pathsFile.split("/")[-1].split(".txt")[0] + ".gfa")

                # Initiate the GFA file.
                if not os.path.exists(gfaFile):
                    with open(gfaFile, "w") as f:
                        gfa = gfapy.Gfa()
                        gfa.add_line("H\tVN:Z:2.0")
                        gfa.to_file(gfaFile)

                # Add corresponding S lines (Segments) to the GFA file, if not already in GFA file.
                with open(gfaFile, "a") as f:
                    gfa = gfapy.Gfa.from_file(gfaFile)
                    if str(scaffoldID) not in gfa.segments:
                        gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(str(scaffoldID), str(len(seq)), os.path.join(outDir, scaffold_fastaFile)))                  
                        gfa.to_file(gfaFile)

            #----------------------------------------------------
            # GFA file - Gap lines
            #----------------------------------------------------
            # Iterate over the links between scaffolds in the PATH.
            for i in range(len(scaffoldsList)-1):

                # Get the orientations of the scaffolds Left and Right
                ## Left
                leftScaffOrient = re.search("[fr]", str(scaffoldsList[i]))
                if leftScaffOrient == None:
                    leftScaffOrient = "fr"
                else:
                    leftScaffOrient = leftScaffOrient.group()
                    if leftScaffOrient == "f":
                        left_orient = "+"
                    elif leftScaffOrient == "r":
                        left_orient = "-"

                ## Right
                rightScaffOrient = re.search("[fr]", str(scaffoldsList[i+1]))
                if rightScaffOrient == None:
                    rightScaffOrient = "fr"
                else:
                    rightScaffOrient = rightScaffOrient.group()
                    if rightScaffOrient == "f":
                        right_orient = "+"
                    elif rightScaffOrient == "r":
                        right_orient = "-"

                ### Particular Case: undetermined orientation (?)
                if (leftScaffOrient == "fr") and (rightScaffOrient == "fr"):

                    # Both left and right scaffolds are forward
                    left_orient = "+"
                    right_orient = "+"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                    # Left scaffold is forward, right scaffold is reverse.
                    left_orient = "+"
                    right_orient = "-"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                    # Left scaffold is reverse, right scaffold is forward.
                    left_orient = "-"
                    right_orient = "+"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                    # Both left and right scaffolds are reverse
                    left_orient = "-"
                    right_orient = "-"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                elif (leftScaffOrient == "fr"):

                    # Left is forward
                    left_orient = "+"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                    # Left is reverse
                    left_orient = "-"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                elif (rightScaffOrient == "fr"):
                    # Right is forward
                    right_orient = "+"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)

                    # Right is reverse
                    right_orient = "-"
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
                        gfa.add_line("G\t*\t{}\t{}\t0\t*".format(str(leftScaffold), str(rightScaffold)))
                        gfa.to_file(gfaFile)
                    
                ### Normal Case: determined orientation (f) or (r)
                else:
                    # Get the names (ID+orient) of the scaffolds Left and Right.
                    leftScaffold = scaffoldsList[i].split("(")[0] + left_orient
                    rightScaffold = scaffoldsList[i+1].split("(")[0] + right_orient

                    #Add corresponding G lines (Gaps) to GFA file (dist=0 means len(gap) unknown)
                    with open(gfaFile, "a") as f:
                        gfa = gfapy.Gfa.from_file(gfaFile)
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

