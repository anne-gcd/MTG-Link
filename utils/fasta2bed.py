#!/usr/bin/env python3
#*****************************************************************************
#  Name: MTG-Link
#  Description: gap-filling tool for draft genome assemblies, dedicated to 
#  linked read data generated by 10XGenomics Chromium technology.
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
import re
import sys
import argparse
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="fasta2bed.py", usage="%(prog)s -in <fasta_file> -out <output_directory>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a FASTA file to a BED file containing the positions of 'Ns' for each scaffold"))

parser.add_argument("-in", dest="input", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-out", dest="outDir", action="store", help="Output directory for saving the BED file", required=True)

args = parser.parse_args()

if (re.match('^.*.fasta$', args.input) is None) and (re.match('^.*.fa$', args.input) is None):
    parser.error("The suffix of the input FASTA file should be: '.fasta' or '.fa'")

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
# FASTA to BED
#----------------------------------------------------
try:
    # Output BED file.
    bedFile = fasta_name.split('.fa')[0] + "_positions_Ns.bed"

    # Iterate over the scaffolds (records) in the FASTA file.
    with open(fasta_file, "r") as in_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            scaffold_name = record.id
            seq = record.seq

            # Get the positions of Ns for each scaffold.
            end = 0
            while "N" in str(seq[end:len(seq)]):
                match = re.search("N+", str(seq[end:len(seq)]))
                unknown_seq = match.group()

                # Positions of the unknown sequence.
                start_unknown_seq = match.start() + end         #0-based, inclusive
                end_unknown_seq = match.end() + end             #non-inclusive

                end = match.end() + end

                # Update the BED file with the positions of the unknown sequence.
                with open(bedFile, "a") as bed_file:
                    bed_file.write("{}\t{}\t{}\n".format(scaffold_name, start_unknown_seq, end_unknown_seq))


except Exception as e:
    print("\nException-")
    print(e)
    sys.exit(1)


print("\nThe BED output file is saved in " + outDir)

