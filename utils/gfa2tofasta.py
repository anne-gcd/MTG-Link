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
parser = argparse.ArgumentParser(prog="gfa2fasta.py", usage="%(prog)s -in <gfaFile.gfa> -out <outDir>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Convert a GFA file (GFA 2.0) to a FASTA file (gaps are returned as 'Ns' regions)"))

parser.add_argument("-in", dest="input", action="store", help="GFA 2.0 file (format: 'xxx.gfa')", required=True)
parser.add_argument("-out", dest="outDir", action="store", help="Directory to save the output FASTA file", required=True)

args = parser.parse_args()

if re.match('^.*.gfa$', args.input) is None:
    parser.error("The suffix of the input GFA file (GFA 2.0) should be: '.gfa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
#GFA 2.0 file
gfa_file = os.path.abspath(args.input)
if not os.path.exists(args.input):
    parser.error("The path of the input GFA file doesn't exist")
gfa_name = gfa_file.split('/')[-1]
print("\nInput GFA file: " + gfa_file,file=sys.stderr )

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
print("The results are saved in " + outDir, file=sys.stderr)

#----------------------------------------------------
# GFA 2.0 to FASTA
#----------------------------------------------------
try:
    segments = {}
    edges = {}
    gaps = []
    paths = []
    added_to_path = []
    fasta_name = gfa_name.split('.gfa')[0] + "_assembly.fasta"

    # Open the input GFA file.
    with open(gfa_file, "r") as f:
        gfa = gfapy.Gfa.from_file(gfa_file)

        # Store the segments in the dict 'segments'.
        for line in gfa.segments:
            if str(line).split('\t')[4] != "co:Z:GFAPY_virtual_line":
                seq_name = str(line).split('\t')[1].split('.k')[0]
                seq_length = str(line).split('\t')[2]

                ## Sequence of the current segment provided directly in the GFA file
                if str(line).split('\t')[3] != "*":
                    seq_sequence = str(line).split('\t')[3]

                    # Append the dict 'segments' with the current segment.
                    segments[seq_name] = [seq_length, seq_sequence]

                ## Sequence of the current segment provided in a FASTA file
                else:
                    seq_file = str(line).split('\t')[4].split(":")[2]
                    with open(seq_file, "r") as input_fasta:
                        for record in SeqIO.parse(input_fasta, "fasta"):
                            if seq_name in str(record.id):
                                seq_sequence = str(record.seq)
                    
                    # Append the dict 'segments' with the current segment.
                    segments[seq_name] = [seq_length, seq_sequence]

        # Store the edges in the dict 'edges'.
        for line in gfa.edges:
            ## s1,s2: seq names
            s1 = str(line).split('\t')[2][:-1]
            s2 = str(line).split('\t')[3][:-1]
            ## o1,o2: orientations of s1,s2 resp.
            o1 = str(line).split('\t')[2][-1:]
            o2 = str(line).split('\t')[3][-1:]
            ## seq1,seq2: seq+orient
            seq1 = str(line).split('\t')[2]
            seq2 = str(line).split('\t')[3]
        
            ## b1: beginning position of s1
            b1 = int(str(line).split('\t')[4])
            ## e1: ending position of s1 (remove the '$' sign if any)
            e1 = str(line).split('\t')[5]
            if (e1[-1:] == '$'):
                e1 = e1[:-1]
            e1 = int(e1) 

            ## b2: beginning position of s2
            b2 = int(str(line).split('\t')[6])
            ## e2: ending position of s2 (remove the '$' sign if any)
            e2 = str(line).split('\t')[7]
            if (e2[-1:] == '$'):
                e2 = e2[:-1]
            e2 = int(e2)

            # Append the dict 'edges' with the current edge.
            if s1 in edges:
                edges[s1].append([s2, o2, b2, e2])
            else:
                edges[s1] = [[o1, b1, e1], [s2, o2, b2, e2]]

        # Store the gaps in the list 'gaps'.
        for line in gfa.gaps:
            ## seq1,seq2: seq+orient
            seq1 = str(line).split('\t')[2]
            seq2 = str(line).split('\t')[3]

            ## g: gap length
            g = str(line).split('\t')[4]
            if g == "*":
                gap_length = "."
            else:
                gap_length = int(g)

            # Append the list 'gaps' with the current gap.
            gaps.append([seq1, seq2, gap_length])


        #################
        # PATHS creation
        #################
        # Iterate over the segments.
        for sname in segments.keys():

            # Current segment has at least one edge and has not been added in a path yet.
            if (sname in edges.keys()) and (sname not in added_to_path):
                ocurr, bcurr, ecurr = [elt for elt in edges[sname][0]]
                
                # Iterate over all the possible edges of the current segment.
                for i in range(1, len(edges[sname])):
                    snext, onext, bnext, enext = [elt for elt in edges[sname][i]]
                    
                    # snext has at least one edge.
                    if snext in edges.keys():
                        ocurrnext, bcurrnext, ecurrnext = [elt for elt in edges[snext][0]]

                        # Iterate over the edge of snext.
                        snextnext, onextnext, bnextnext, enextnext = [elt for elt in edges[snext][1]]

                        # Path creation.
                        overlap1 = ecurr - bcurr
                        overlap2 = ecurrnext - bcurrnext
                        path = str(sname)+str(ocurr) +","+ str(snext)+str(onext) +","+ str(snextnext)+str(onextnext) +"\t"+ str(overlap1) +"M,"+ str(overlap2) +"M"
                        added_to_path.append(sname)
                        added_to_path.append(snext)
                        added_to_path.append(snextnext)

                        # Append the list 'paths' with the current path.
                        paths.append(path)

            # Current segment has no edge (gap).
            else:
                for gap in gaps:
                    if (gap[0] == str(sname)+"+") or (gap[0] == str(sname)+"-"):
                        seq1 = gap[0]
                        seq2 = gap[1]
                        gap_length = gap[2]
                        if gap_length == ".":
                            g = "."
                        else:
                            g = str(gap_length) + "J"
                        path = str(seq1) +";"+ str(seq2) +"\t"+ str(g)

                        # Append the list 'paths' with the current path.
                        paths.append(path)


        #################
        # FASTA creation
        #################
        # Iterate over the paths.
        for path in paths:

            ## Assemblies
            if "," in path:
                scaffolds = str(path).split('\t')[0].split(',')
                overlaps = str(path).split('\t')[1].split(',')
                assembly = ""

                # Iterate over the different scaffolds of the path.
                for i in range(len(scaffolds)):

                    strand = scaffolds[i][-1:]
                    ## fwd strand
                    if strand == "+":
                        sequence = segments[scaffolds[i][:-1]][1]

                    ## rev strand
                    if strand == "-":
                        sequence = rc(segments[scaffolds[i][:-1]][1])
                    
                    ##Attn: remove overlap only on scaffolds, not on gapfilled seq
                    # Initiation.
                    if i == 0:
                        overlap = int(overlaps[i][:-1])
                        assembly += sequence[:-overlap]

                    # Last scaffold.
                    elif i == len(scaffolds)-1:
                        overlap = int(overlaps[i-1][:-1])
                        assembly += sequence[overlap:]
                    
                    # Assembled sequence.
                    elif i % 2 == 1:
                        assembly += sequence

                    # Scaffold.
                    else:
                        overlap = int(overlaps[i][:-1])
                        assembly += (sequence[overlap:])[:-overlap]

                # Write the assembly sequence to the FASTA file.
                name = ""
                for i in range(1, len(scaffolds)):
                    name += scaffolds[i] + ":"
                name = name[:-1]    #to remove the last ':' char
                name += " len_" + str(len(assembly))

                with open(fasta_name, "a") as fasta:
                    fasta.write(">" + name)
                    fasta.write("\n" + str(assembly) + "\n")

            ## Gaps
            elif ";" in path:
                scaffolds = str(path).split('\t')[0].split(';')
                overlap = str(path).split('\t')[1]

                # Iterate over the different scaffolds of the path.
                for i in range(len(scaffolds)):

                    strand = scaffolds[i][-1:]
                    ## fwd strand
                    if strand == "+":
                        sequence = segments[scaffolds[i][:-1]][1]

                    ## rev strand
                    if strand == "-":
                        sequence = rc(segments[scaffolds[i][:-1]][1])

                    # Write the unassembled sequence to the FASTA file.
                    name = scaffolds[i][:-1]
                    name += " len_" + str(len(sequence))

                    with open(fasta_name, "a") as fasta:
                        fasta.write(">" + name)
                        fasta.write("\n" + str(sequence) + "\n")


except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe output FASTA file is saved in " + outDir)
print("Output FASTA file: " + fasta_name)

