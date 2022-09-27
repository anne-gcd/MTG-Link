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
import gfapy
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="bed2gfa.py", usage="%(prog)s -bed <bedFile.bed> -fa <fastaFile.fasta> -out <outputGFAFile> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a BED file containing the 'N's coordinates for each scaffold (or locus coordinates) to a GFA file (GFA 2.0) ('N's regions are treated as gaps). We can filter the 'N's regions by their size (e.g. gap lengths) and by the contigs' sizes on both sides (long enough for ex to get enough barcodes)"))

parser.add_argument("-bed", dest="bed", action="store", help="BED file containing the 'Ns' coordinates for each scaffold (format: 'xxx.bed')", required=True)
parser.add_argument("-fa", dest="fasta", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-min", dest="min", action="store", type=int, help="Minimum size of the 'Ns' region to treat as a gap")
parser.add_argument("-max", dest="max", action="store", type=int, help="Maximum size of the 'Ns' region to treat as a gap")
parser.add_argument("-contigs", dest="contigs_size", action="store", type=int, help="Minimum size of the flanking contigs of the 'Ns' region to treat as a gap")
parser.add_argument("-out", dest="outGFA", action="store", help="Name of the output GFA file")

args = parser.parse_args()

if re.match('^.*.bed$', args.bed) is None:
    parser.error("The suffix of the input BED file should be: '.bed'")

if (re.match('^.*.fasta$', args.fasta) is None) and (re.match('^.*.fa$', args.fasta) is None):
    parser.error("The suffix of the input FASTA file should be: '.fasta' or '.fa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
# BED file.
bedFile = os.path.abspath(args.bed)
if not os.path.exists(args.bed):
    parser.error("The path of the input BED file doesn't exist")
bed_name = bedFile.split('/')[-1]
print("\nInput BED file: " + bedFile)

# FASTA file.
fastaFile = os.path.abspath(args.fasta)
if not os.path.exists(args.fasta):
    parser.error("The path of the input FASTA file doesn't exist")
fasta_name = fastaFile.split('/')[-1]
print("Input FASTA file: " + fastaFile)

#----------------------------------------------------
# Directory for saving results
#----------------------------------------------------
cwd = os.getcwd()
outDir = cwd
print("\nThe results are saved in " + outDir)

#----------------------------------------------------
# BED to GFA
#----------------------------------------------------
try:
    # Create a dictionary containing the positions of 'Ns' for each scaffold.
    positions_NsDict = {}

    # Iterate over the positions of 'Ns' in the BED file.
    with open(bedFile, "r") as bed_file:
        for line in bed_file:
            chrom = line.split('\t')[0]
            chromStart = line.split('\t')[1] 
            chromEnd = line.split('\t')[2].split('\n')[0]

            # Get the size of the 'Ns' region (e.g. the gap size).
            gap_size = int(chromEnd) - int(chromStart)

            # If 'min' < gap_size < 'max', append the dictionary 'positions_NsDict' with the 'Ns' regions (e.g. gaps) for each scaffold.
            ## key = scaffold_name ; value = list of positions of gaps
            if (args.min is not None) and (gap_size < int(args.min)):
                continue
            if (args.max is not None) and (gap_size > int(args.max)):
                continue
            if chrom in positions_NsDict:
                positions_NsDict[chrom].append([int(chromStart), int(chromEnd)])
            else:
                positions_NsDict[chrom] = [[int(chromStart), int(chromEnd)]]
        
    # For each gap, get the gap's sequence and the left and right flanking sequences.
    with open(fastaFile, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in positions_NsDict:

                # Iterate over the gaps of 'positions_NsDict' for each scaffold.
                gap_coordinatesList = positions_NsDict[record.id]
                gap_coordinatesList.sort()
                for i in range(len(gap_coordinatesList)):

                    # Get the length of the gap.
                    gap_length = int(gap_coordinatesList[i][1]) - int(gap_coordinatesList[i][0])

                    # Get the left flanking sequence.
                    if i == 0:      #first gap
                        left_start_index = 0 
                        left_end_index = gap_coordinatesList[i][0]
                    else:
                        left_start_index = gap_coordinatesList[i-1][1]
                        left_end_index = gap_coordinatesList[i][0]
                    left_flanking_seq = record.seq[int(left_start_index):int(left_end_index)]

                    # Get the right flanking sequence.
                    if i == len(gap_coordinatesList)-1:     #last gap
                        right_start_index = gap_coordinatesList[i][1]
                        right_end_index = len(record.seq)
                    else:
                        right_start_index = gap_coordinatesList[i][1]
                        right_end_index = gap_coordinatesList[i+1][0]
                    right_flanking_seq = record.seq[int(right_start_index):int(right_end_index)]

                    # If left and right sequences > 'contigs_size', add the gap to the output GFA file.
                    if (args.contigs_size is not None) and (len(left_flanking_seq) < int(args.contigs_size)) and (len(right_flanking_seq) < int(args.contigs_size)):
                        continue
                    else:
                        os.chdir(outDir)

                        #----------------------------------------------------
                        # FASTA files (left/right sequences)
                        #----------------------------------------------------
                        if args.min is None:
                            minLength = "noMin"
                        else:
                            minLength = args.min
                        if args.max is None:
                            maxLength = "noMax"
                        else:
                            maxLength = args.max
                        if args.contigs_size is None:
                            minContigSize = "noMinFlank"
                        else:
                            minContigSize = args.contigs_size

                        # Left flanking contig.
                        left_name = str(record.id) +"_"+ str(left_start_index) +"-"+ str(left_end_index) +"-L"
                        left_contig = left_name + "+"
                        left_length = len(left_flanking_seq)

                        # Right flanking contig.
                        right_name = str(record.id) +"_"+ str(right_start_index) +"-"+ str(right_end_index) +"-R"
                        right_contig = right_name + "+"
                        right_length = len(right_flanking_seq)

                        # Output FASTA files containing the flanking sequences of the current gap.
                        ##Left
                        left_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_" + str(record.id) +"_"+ str(left_start_index) +"-"+ str(left_end_index) +".g"+ str(gap_length) + ".left.fasta")
                        with open(left_fastaFile, "w") as left:
                            left.write(">{} _ len {}".format(left_name, str(left_length)))
                            left.write("\n" + str(left_flanking_seq) + "\n")

                        ##Right
                        right_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_" + str(record.id) +"_"+ str(right_start_index) +"-"+ str(right_end_index) +".g"+ str(gap_length) + ".right.fasta")
                        with open(right_fastaFile, "w") as right:
                            right.write((">{} _ len {}".format(right_name, str(right_length))))
                            right.write("\n" + str(right_flanking_seq) + "\n")
                        
                        #----------------------------------------------------
                        # GFA file
                        #----------------------------------------------------
                        os.chdir(outDir)
                        if args.outGFA is not None:
                            gfaFile = os.path.abspath(str(args.outGFA))
                        else:
                            gfaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_gaps_" + str(minLength) + "-" + str(maxLength) + "_contigs_" + str(minContigSize) + ".gfa")

                        # Initiate the GFA file.
                        if not os.path.exists(gfaFile):
                            with open(gfaFile, "w") as f:
                                gfa = gfapy.Gfa()
                                gfa.add_line("H\tVN:Z:2.0")
                                gfa.to_file(gfaFile)

                        # Add corresponding lines to the GFA file.
                        with open(gfaFile, "a") as f:
                            gfa = gfapy.Gfa.from_file(gfaFile)
                            gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(left_name, left_length, os.path.join(outDir, left_fastaFile)))
                            gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(right_name, right_length, os.path.join(outDir, right_fastaFile)))
                            gfa.add_line("G\t*\t{}\t{}\t{}\t*".format(left_contig, right_contig, gap_length))
                            gfa.to_file(gfaFile)


except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe output GFA file and the corresponding FASTA files are saved in " + outDir)
print("Output GFA file: " + gfaFile)

