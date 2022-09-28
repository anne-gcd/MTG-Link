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
import re
import sys
import argparse
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="fasta2bed.py", usage="%(prog)s -fa <fastaFile.fasta> -out <outputBEDFile>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a FASTA file to a BED file containing the 'N's coordinates for each scaffold"))

parser.add_argument("-fa", dest="fasta", action="store", help="FASTA file containing the sequences of the scaffolds (reference genome) (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-out", dest="outBED", action="store", help="Name of the output BED file")

args = parser.parse_args()

if (re.match('^.*.fasta$', args.fasta) is None) and (re.match('^.*.fa$', args.fasta) is None):
    parser.error("The suffix of the input FASTA file should be: '.fasta' or '.fa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
fastaFile = os.path.abspath(args.fasta)
if not os.path.exists(args.fasta):
    parser.error("The path of the input FASTA file doesn't exist")
fasta_name = fastaFile.split('/')[-1]
print("\nInput FASTA file: " + fastaFile)

#----------------------------------------------------
# Directory for saving results
#----------------------------------------------------
cwd = os.getcwd()
outDir = cwd
print("\nThe results are saved in " + outDir)

#----------------------------------------------------
# FASTA to BED
#----------------------------------------------------
try:
    # Output BED file.
    if args.outBED is not None:
        bedFile = str(args.outBED)
    else:
        bedFile = fasta_name.split('.fa')[0] + "_NsCoord.bed"

    # Iterate over the scaffolds (records) in the FASTA file.
    with open(fastaFile, "r") as in_fasta:
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
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe BED output file is saved in " + outDir)
print("Output BED file: " + bedFile)

