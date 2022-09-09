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


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="mergegfa.py", usage="%(prog)s -1 <gfaFile1.gfa> -2 <gfaFile2.gfa> -out <merged_gfa_file>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Merge two GFA files together"))

parser.add_argument("-1", dest="gfa1", action="store", help="GFA 2.0 file n°1 (format: 'xxx.gfa')", required=True)
parser.add_argument("-2", dest="gfa2", action="store", help="GFA 2.0 file n°2 (format: 'xxx.gfa')", required=True)
parser.add_argument("-out", dest="merged_gfa", action="store", help="Name of the output merged GFA file", required=True)

args = parser.parse_args()

if re.match('^.*.gfa$', args.gfa1) is None:
    parser.error("The suffix of the input GFA file n°1 (GFA 2.0) should be: '.gfa'")
if re.match('^.*.gfa$', args.gfa2) is None:
    parser.error("The suffix of the input GFA file n°2 (GFA 2.0) should be: '.gfa'")
#----------------------------------------------------
# Input files
#----------------------------------------------------
#GFA 2.0 file no.1 
gfa_file1 = os.path.abspath(args.gfa1)
if not os.path.exists(gfa_file1):
    parser.error("The path of the input GFA file n°1 (GFA 2.0)" + args.gfa1 + "doesn't exist")

#GFA 2.0 file no.2
gfa_file2 = os.path.abspath(args.gfa2)
if not os.path.exists(gfa_file2):
    parser.error("The path of the input GFA file n°2 (GFA 2.0)" + args.gfa2 + "doesn't exist")

#----------------------------------------------------
# Merge the two GFA files together
#----------------------------------------------------
try:
    S_seen = {}
    E_seen = {}
    G_seen = {}

    #Open the two input GFA files
    gfa1 = gfapy.Gfa.from_file(gfa_file1)
    gfa2 = gfapy.Gfa.from_file(gfa_file2)

    #Output merged GFA file
    merged_gfa_file = str(args.merged_gfa)
    with open(merged_gfa_file, "w") as f:
        merged_gfa = gfapy.Gfa()
        merged_gfa.add_line("H\tVN:Z:2.0")

        #Iterate over the 'S' lines of the input GFA file no.1 and add them to the output merged GFA
        for line in gfa1.segments:
            seq_name = str(line).split('\t')[1]
            S_seen[seq_name]=1
            merged_gfa.add_line(str(line))

        #Iterate over the 'S' lines of the input GFA file no.2
        for line in gfa2.segments:
            seq_name = str(line).split('\t')[1]
            #Add the non already present 'S' lines to the output merged GFA
            if not seq_name in S_seen:
                merged_gfa.add_line(str(line))

        #Iterate over the 'E' lines of the input GFA file no.1 and add them to the output merged GFA
        for line in gfa1.edges:
            edge_name = str(line).split('\t')[2] + "--" + str(line).split('\t')[3]
            E_seen[edge_name]=1
            merged_gfa.add_line(str(line))

        #Iterate over the 'E' lines of the input GFA file no.2
        for line in gfa2.edges:
            edge_name = str(line).split('\t')[2] + "--" + str(line).split('\t')[3]
            #Add the non already present 'E' lines to the output merged GFA
            if not edge_name in E_seen:
                merged_gfa.add_line(str(line))

        #Iterate over the 'G' lines of the input GFA file no.1 and add them to the output merged GFA
        for line in gfa1.gaps:
            gap_name = str(line).split('\t')[2] + "--" + str(line).split('\t')[3]
            G_seen[gap_name]=1
            merged_gfa.add_line(str(line))
        
        #Iterate over the 'G' lines of the input GFA file no.2
        for line in gfa2.gaps:
            gap_name = str(line).split('\t')[2] + "--" + str(line).split('\t')[3]
            if not gap_name in G_seen:
                merged_gfa.add_line(str(line))

        merged_gfa.to_file(merged_gfa_file)


except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print("\nException-")
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("Output merged GFA file: " + merged_gfa_file)

