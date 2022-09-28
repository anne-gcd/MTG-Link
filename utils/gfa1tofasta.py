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
parser = argparse.ArgumentParser(prog="gfa1tofasta.py", usage="%(prog)s -in <gfa1File.gfa> -out <outDir>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a GFA file (GFA 1.0) to a FASTA file (gaps are returned as 'Ns' regions)"))

parser.add_argument("-in", dest="input", action="store", help="GFA 1.0 file (format: 'xxx.gfa')", required=True)
parser.add_argument("-out", dest="outDir", action="store", help="Directory to save the output FASTA file", required=True)

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
    i = 0
    fasta_name = gfa_name.split('.gfa')[0] + "_assembly.fasta"

    with open(gfa_file, "r") as f:
        gfa = gfapy.Gfa.from_file(gfa_file)

        # Iterate over the 'Segment' lines of the input GFA.
        for line in gfa.segments:
            print(line)
            if str(line).split('\t')[4] != "co:Z:GFAPY_virtual_line":
                seq_name = str(line).split('\t')[1]
                seq_file = str(line).split('\t')[4].split(":")[2]

                # If the segment's sequence is not in the list 'fasta_links', add it to this list.
                ## to open the specified 'seq_file' only once (even if it appears for several segments in GFA).
                if seq_file not in fasta_links:
                    fasta_links.append(seq_file)

        #Dictionary containing the scaffolds sequences and the gapfilled sequences
        for seq_file in fasta_links:
            if seq_file == "GFAPY_virtual_line":
                continue
            else:
                with open(seq_file, "r") as input_fasta:
                    for record in SeqIO.parse(input_fasta, "fasta"):
                        name = (record.id).split(" ")[0]
                        seq = str(record.seq)
                        all_seq[name] = seq
                  
        #Get the paths from the 'Path' line of the input GFA
        #(only the fwd gapfilled seq were kept to construct the path)
        for line in gfa.paths:
            line = str(line).split('\t')
            scaffolds = str(line[2]).split(',')
            overlap_lengths = str(line[3]).split(',')

            #Iterate over the different scaffolds of the path
            assembly = ""
            for i in range(len(scaffolds)):

                #fwd strand
                strand = scaffolds[i][-1:]
                if strand == "+":
                    sequence = all_seq[scaffolds[i][:-1]]
                #rev strand
                if strand == "-":
                    sequence = rc(all_seq[scaffolds[i][:-1]])

                #Attn: remove overlap only on scaffolds, not on gapfilled seq
                #initiation
                if i == 0:
                    over = int(overlap_lengths[i][:-1])
                    assembly += sequence[:-over]
                
                #last scaffold
                elif i == len(scaffolds)-1:
                    over = int(overlap_lengths[i-1][:-1])
                    assembly += sequence[over:]

                #if gapfilled sequence
                elif i % 2 == 1:
                    assembly += sequence

                #if scaffold
                else:
                    over = int(overlap_lengths[i][:-1])
                    assembly += (sequence[over:])[:-over]

            #Write the assembly sequence to the FASTA file
            name = ""
            for i in range(0, len(scaffolds), 2):
                name += scaffolds[i][:-1]
                if scaffolds[i][-1:] == "+":
                    name += "f_"
                else:
                    name += "r_"
            name = name[:-1]    #to remove the last '_' char
            name += " len " + str(len(assembly))

            with open(fasta_name, "a") as fasta:
                fasta.write(">" + name)
                fasta.write("\n" + str(assembly) + "\n")
            

except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe FASTA output file is saved in " + outDir)


#TODO: return gaps as Ns regions ?? 