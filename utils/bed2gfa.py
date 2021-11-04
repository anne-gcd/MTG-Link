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
import gfapy
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="bed2fasta.py", usage="%(prog)s -bed <bed_file> -fa <fasta_file> -out <output_directory> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a BED file containing the positions of 'Ns' for each scaffold to a GFA file ('Ns' regions are treated as gaps). We can filter the 'Ns' regions by their size (e.g. gap sizes) and by the contigs' sizes on both sides (long enough for ex to get enough barcodes)"))

parser.add_argument("-bed", dest="bed", action="store", help="BED file containing the positions of 'Ns' for each scaffold (format: 'xxx.bed')", required=True)
parser.add_argument("-fa", dest="fasta", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-min", dest="min", action="store", type=int, help="Minimum size of the 'Ns' region to treat/process as a gap")
parser.add_argument("-max", dest="max", action="store", type=int, help="Maximum size of the 'Ns' region to treat/process as a gap")
parser.add_argument("-contigs", dest="contigs_size", action="store", type=int, help="Minimum size of the flanking contigs of the 'Ns' region to treat/process as a gap")
parser.add_argument("-out", dest="outDir", action="store", help="Output directory for saving the GFA file", required=True)

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
print("\nInput FASTA file: " + fastaFile)

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
            chromEnd = line.split('\t')[2]

            # Get the size of the 'Ns' region (e.g. the gap size).
            gap_size = int(chromEnd) - int(chromStart)

            # If 'min' < gap_size < 'max', append the dictionary 'positions_NsDict' with the 'Ns' regions (e.g. gaps) for each scaffold.
            ## key = scaffold_name ; value = list of positions of gaps
            if (args.min is not None) and (gap_size < args.min):
                break
            if (args.max is not None) and (gap_size > args.max):
                break
            if chrom in positions_NsDict:
                positions_NsDict[chrom].append([chromStart, chromEnd])
            else:
                positions_NsDict[chrom] = [[chromStart, chromEnd]]
        
    # For each gap, get the gap's sequence and the left and right flanking sequences.
    with open(fastaFile, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in positions_NsDict:

                # Iterate over the gaps of 'positions_NsDict' for each scaffold:
                gap_coordinatesList = positions_NsDict[record.id]
                gap_count = 0
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


                    # If left and right sequences > 'contigs_size', continue (e.g. add the gap to the output GFA file).
                    if (args.contigs_size is not None) and (len(left_flanking_seq) < args.contigs_size) and (len(right_flanking_seq) < args.contigs_size):
                        break
                    else:
                        gap_count += 1
                        os.chdir(outDir)

                        #----------------------------------------------------
                        # FASTA files (left/right sequences)
                        #----------------------------------------------------
                        if args.min is None:
                            args.min = "."
                        if args.max is None:
                            args.max = "."
                        if args.contigs_size is None:
                            args.contigs_size = "."

                        # Output FASTA file containing the flanking sequences of the gap regions of the output GFA file.
                        out_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_flanking_seq_gaps_" + str(args.min) + "-" + str(args.max) + "_contigs_" + str(args.contigs_size) + ".fasta")
                        print("\nOutput FASTA file: " + out_fastaFile)
                                    
                        # Save the left and right sequences to the output FASTA file.
                        left_name = record.id + "_gap" + str(gap_count) + "-L"
                        left_contig = left_name + "+"
                        left_length = len(left_flanking_seq)  

                        right_name = record.id + "_gap" + str(gap_count) + "-R"
                        right_contig = right_name + "+"
                        right_length = len(right_flanking_seq)

                        with open(out_fastaFile, "a") as output:
                            output.write(">{} _ len {}".format(left_name, left_length))
                            output.write("\n" + str(left_flanking_seq) + "\n")
                            output.write(">{} _ len {}".format(right_name, right_length))
                            output.write("\n" + str(right_flanking_seq) + "\n")

                        #----------------------------------------------------
                        # GFA file
                        #----------------------------------------------------
                        os.chdir(outDir)

                        gfaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_gaps_" + str(args.min) + "-" + str(args.max) + "_contigs_" + str(args.contigs_size) + ".gfa")
                        print("\nOutput GFA file: " + gfaFile)

                        # Initiate the GFA file.
                        if not os.path.exists(gfaFile):
                            with open(gfaFile, "w") as f:
                                gfa = gfapy.Gfa()
                                gfa.add_line("H\tVN:Z:2.0")
                                gfa.to_file(gfaFile)

                        # Add corresponding lines to the GFA file.
                        with open(gfaFile, "a") as f:
                            gfa = gfapy.Gfa.from_file(gfaFile)
                            gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(left_name, left_end_index, os.path.join(outDir, out_fastaFile)))
                            gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(right_name, right_start_index, os.path.join(outDir, out_fastaFile)))
                            gfa.add_line("G\t*\t{}\t{}\t{}\t*".format(left_contig, right_contig, gap_length))
                            gfa.to_file(gfaFile)


except Exception as e:
    print("\nException-")
    print(e)
    sys.exit(1)


print("\nThe output GFA file and the corresponding FASTA files are saved in " + outDir)

