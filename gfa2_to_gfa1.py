#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import re
import argparse
import gfapy


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="gfa.2_to_gfa.1.py", usage="%(prog)s -in <input_gfa_2.0)> -ext <extension_gap_size> -out <output_directory>", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Convert a GFA 2.0 file into a GFA 1.0 file
                                '''))

parser.add_argument("-in", "--input", action="store", help="GFA 2.0 file (format: 'xxx.gfa')", required=True)
parser.add_argument("-ext", "--extension", action="store", type=int, help="size of the extension of the gap", required=True)
parser.add_argument("-out", "--outDir", action="store", help="output directory for saving the GFA 1.0 file", required=True)

args = parser.parse_args()

if re.match('^.*.gfa$', args.input) is None:
    parser.error("The suffix of the input GFA file should be: '.gfa'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
input_gfa = os.path.abspath(args.input)
if not os.path.exists(args.input):
    parser.error("The path of the input GFA file doesn't exist")
gfa_name = (input_gfa.split('/')[-1]).split('.gfa')[0]
print("\nInput GFA file: " + input_gfa)

ext = args.extension

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
# GFA 2.0 to GFA 1.0
#----------------------------------------------------
try:
    output_gfa = outDir + "/" + gfa_name + "_1.0.gfa"

    #Iterate over the 'Edge' lines of the input GFA (GFA 2.0)
    with open(input_gfa, "r") as f2, open(output_gfa, "a") as f1:
        gfa2 = gfapy.Gfa.from_file(input_gfa)
        gfa1 = gfapy.Gfa.from_file(output_gfa)

        gfa1.add_line("H\tVN:Z:1.0")

        for line in gfa2.segments:

            name = str(line).split('\t')[1]
            length = str(line).split('\t')[2]
            sequence_file = str(line).split('\t')[4]

            gfa1.add_line("S\t{}\t*\tLN:i:{}\t{}".format(name, length, sequence_file))

        for line in gfa2.edges:

            s1_orient = str(line).split('\t')[2]
            s1 = "".join(list(s1_orient)[:-1])
            orient1 = list(s1_orient)[-1]
            
            s2_orient = str(line).split('\t')[3]
            s2 = "".join(list(s2_orient)[:-1])
            orient2 = list(s2_orient)[-1]

            overlap_length = int((str(line).split('\t')[-2]).split('$')[0]) - int(str(line).split('\t')[-3]) + 1

            gfa1.add_line("L\t{}\t{}\t{}\t{}\t{}M".format(s1, orient1, s2, orient2, overlap_length))

        gfa1.to_file(output_gfa)


except Exception as e:
    print("\nException-")
    exc_type, exc_tb = sys.exc_info()
    print(exc_type, exc_tb.tb_lineno)
    sys.exit(1)