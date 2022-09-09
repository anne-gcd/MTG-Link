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
parser = argparse.ArgumentParser(prog="vcf2gfa.py", usage="%(prog)s -vcf <vcfFile.vcf> -fa <fastaFile.fasta> -out <outDir> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a VCF file containing the insertions coordinates to a GFA file (GFA 2.0). We can filter the insertions by the contigs' sizes on both sides (long enough for ex to get enough barcodes). As for some insertions, there are a few micro-homologies (mh), we are not sure that the given position in the VCF file is the exact position, so extend the insertion by '--extension' bp on both sides of the insertion"))

parser.add_argument("-vcf", dest="vcf", action="store", help="VCF file containing the insertions coordinates (format: 'xxx.vcf')", required=True)
parser.add_argument("-fa", dest="fasta", action="store", help="FASTA file containing the sequences of the scaffolds obtained from the assembly (format: 'xxx.fasta' or 'xxx.fa')", required=True)
parser.add_argument("-contigs", dest="contigs_size", action="store", type=int, help="Minimum size of the flanking contigs of the insertion to treat as a target")
parser.add_argument("-extension", dest="ext_size", action="store", type=int, default=50, help="Size of the extension of the insertion (bp): extend the insertion by '--extension' bp on both sides of the insertion [default: 50]")
parser.add_argument("-out", dest="outDir", action="store", help="Directory to save the output GFA file and gap flanking sequences FASTA files", required=True)

args = parser.parse_args()

if re.match('^.*.vcf$', args.vcf) is None:
    parser.error("The suffix of the input VCF file should be: '.vcf'")

if (re.match('^.*.fasta$', args.fasta) is None) and (re.match('^.*.fa$', args.fasta) is None):
    parser.error("The suffix of the input FASTA file should be: '.fasta' or '.fa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
# VCF file.
vcfFile = os.path.abspath(args.vcf)
if not os.path.exists(args.vcf):
    parser.error("The path of the input VCF file doesn't exist")
vcf_name = vcfFile.split('/')[-1]
print("\nInput VCF file: " + vcfFile)

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
# VCF to GFA
#----------------------------------------------------
try:
    # Create a dictionary containing the insertions' coordinates for each scaffold.
    ins_Coords = {}

    # Iterate over the positions of 'Ns' in the BED file.
    with open(vcfFile, "r") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            else:
                chrom, pos, sv_id, ref, alt, __, __, info, *__ = line.rstrip().split("\t")
                sv_type = str(info).split('SVTYPE=')[1].split(';')[0]
                
                # Process only the insertions (SVTYPE=INS).
                if sv_type == "INS":

                    ## Get the start and end positions of the insertion' coordinates.
                    insCoord_start = int(pos) - args.ext_size
                    end_pos = str(info).split('END=')[1]
                    if ";" in end_pos:
                        end_pos = end_pos.split(';')[0]
                    insCoord_end = int(end_pos) + args.ext_size
                    
                    ## Get the size of the insertion region (e.g. the target size).
                    target_size = insCoord_end - insCoord_start

                    ## Append the dictionary 'ins_Coords' with the insertions' coordinates (e.g. targets) for each scaffold.
                    ## key = scaffold_name ; value = list of positions of targets
                    if chrom in ins_Coords:
                        ins_Coords[chrom].append([insCoord_start, insCoord_end])
                    else:
                        ins_Coords[chrom] = [[insCoord_start, insCoord_end]]

    # For each insertion, get the insertion's length and the left and right flanking sequences.
    with open(fastaFile, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in ins_Coords:

                # Iterate over the insertions of 'ins_Coords' for each scaffold.
                insertions_coordinatesList = ins_Coords[record.id]
                insertions_coordinatesList.sort()
                for i in range(len(insertions_coordinatesList)):

                    # Get the length of the insertion.
                    insertion_length = int(insertions_coordinatesList[i][1]) - int(insertions_coordinatesList[i][0])

                    # Get the left flanking sequence.
                    if i == 0:      #first target
                        left_start_index = 0 
                        left_end_index = insertions_coordinatesList[i][0]
                    else:
                        left_start_index = insertions_coordinatesList[i-1][1]
                        left_end_index = insertions_coordinatesList[i][0]
                    left_flanking_seq = record.seq[int(left_start_index):int(left_end_index)]

                    # Get the right flanking sequence.
                    if i == len(insertions_coordinatesList)-1:     #last target
                        right_start_index = insertions_coordinatesList[i][1]
                        right_end_index = len(record.seq)
                    else:
                        right_start_index = insertions_coordinatesList[i][1]
                        right_end_index = insertions_coordinatesList[i+1][0]
                    right_flanking_seq = record.seq[int(right_start_index):int(right_end_index)]

                    # If left and right sequences > 'contigs_size', add the insertion/target to the output GFA file.
                    if (args.contigs_size is not None) and (len(left_flanking_seq) < int(args.contigs_size)) and (len(right_flanking_seq) < int(args.contigs_size)):
                        continue
                    else:
                        os.chdir(outDir)

                        #----------------------------------------------------
                        # FASTA files (left/right sequences)
                        #----------------------------------------------------
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
                        left_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_" + str(record.id) +"_"+ str(left_start_index) +"-"+ str(left_end_index) +".g"+ str(insertion_length) +".c"+ str(minContigSize) + ".left.fasta")
                        with open(left_fastaFile, "w") as left:
                            left.write(">{} _ len {}".format(left_name, str(left_length)))
                            left.write("\n" + str(left_flanking_seq) + "\n")

                        ##Right
                        right_fastaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_" + str(record.id) +"_"+ str(right_start_index) +"-"+ str(right_end_index) +".g"+ str(insertion_length) +".c"+ str(minContigSize) + ".right.fasta")
                        with open(right_fastaFile, "w") as right:
                            right.write((">{} _ len {}".format(right_name, str(right_length))))
                            right.write("\n" + str(right_flanking_seq) + "\n")
                        
                        #----------------------------------------------------
                        # GFA file
                        #----------------------------------------------------
                        os.chdir(outDir)
                        gfaFile = os.path.abspath(fasta_name.split(".fa")[0] + "_insertions_extension_" + str(args.ext_size) + "_contigs_" + str(minContigSize) + ".gfa")

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
                            gfa.add_line("G\t*\t{}\t{}\t{}\t*".format(left_contig, right_contig, insertion_length))
                            gfa.to_file(gfaFile)


except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe output GFA file and the corresponding FASTA files are saved in " + outDir)
print("Output GFA file: " + gfaFile)

