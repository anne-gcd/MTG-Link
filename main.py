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

"""Module 'main.py': Initialization

The module 'main.py' enables to get the input parameters and variables as well as the directories in which to save the results.
"""

from __future__ import print_function
import argparse
import os
import re
import sys
from Bio import SeqIO


#----------------------------------------------------
# Arg parser
#----------------------------------------------------

MTGLINK_VERSION = "2.4.0"

parser = argparse.ArgumentParser(prog="mtglink.py", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Local assembly with linked read data, using either a De Bruijn Graph (DBG) algorithm or an Iterative Read Overlap (IRO) algorithm"))

parser.add_argument('-v', action="version", version='%(prog)s {version}'.format(version=MTGLINK_VERSION))

subparsers = parser.add_subparsers(dest="MTGLink_module", help="MTGLink module used for the Local Assembly step", required=True)

# DBG module options.
parserDBG = subparsers.add_parser("DBG", help="Local assembly using a De Bruijn Graph (DBG) algorithm")
parserDBG.add_argument('-gfa', dest="gfa", action="store", help="Input GFA file (GFA 2.0) (format: xxx.gfa)", required=True)
parserDBG.add_argument('-bam', dest="bam", action="store", help="BAM file: linked-reads mapped on reference genome (format: xxx.bam)", required=True)
parserDBG.add_argument('-fastq', dest="reads", action="store", help="File (gzipped) of indexed reads (format: xxx.fastq.gz | xxx.fq.gz)", required=True)
parserDBG.add_argument('-index', dest="index", action="store", help="Barcodes index file (format: xxx.bci)", required=True)
parserDBG.add_argument('-out', dest="outDir", action="store", default="./MTG-Link_results", help="Output directory [default: ./MTG-Link_results]")
parserDBG.add_argument('-line', dest="line", action="store", type=int, help="Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]")
parserDBG.add_argument('-bxuDir', dest="bxuDir", action="store", help="Directory where the FASTQ files containing the subsample of reads are located (1 file per target) (format of FASTQ files: xxx.bxu.fastq) [to provide if the read subsampling step has already been done for this dataset]")
parserDBG.add_argument('-t', dest="threads", action="store", type=int, default=1, help="Number of threads to use for the Read Subsampling step [default: 1]")
parserDBG.add_argument('-flank', dest="flankSize", action="store", type=int, default=10000, help="Flanking sequences' size (bp) [default: 10000]")
parserDBG.add_argument('-occ', dest="minBarcOcc", action="store", type=int, default=2, help="Minimum number of occurrences in target flanking regions for a barcode to be retained in the union set [default: 2]")
parserDBG.add_argument('-ext', dest="extSize", action="store", type=int, default=500, help="Size of the extension of the target on both sides (bp); determine start/end of local assembly [default: 500]")
parserDBG.add_argument('-l', dest="maxLength", action="store", type=int, default=10000, help="Maximum assembly length (bp) (it could be a bit bigger than the length of the target to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]")
parserDBG.add_argument('-m', dest="minLength", action="store", type=int, default=1000, help="Minimum assembly length (bp), by default 2*(-ext) bp [default: 1000]")
parserDBG.add_argument('-k', dest="kmerSize", action="store", type=int, default=[51, 41, 31, 21],  nargs='+', help="K-mer size(s) used for the local assembly [default: [51, 41, 31, 21]]")
parserDBG.add_argument('-a', dest="abundanceThreshold", action="store", type=int, default=[3, 2], nargs='+', help="Minimal abundance threshold for solid k-mers [default: [3, 2]]")
parserDBG.add_argument("--force", action="store_true", help="To force search on all '-k' values provided")
parserDBG.add_argument('-max-nodes', dest="maxNodes", action="store", type=int, default=1000, help="Maximum number of nodes in contig graph [default: 1000]")
parserDBG.add_argument('-nb-cores', dest="nbCores", action="store", type=int, default=1, help="Number of cores for the Local Assembly step (DBG assembly) [default: 1]")
parserDBG.add_argument('-max-memory', dest="maxMemory", action="store", type=int, default=0, help="Maximum memory for graph building (in MBytes) [default: 0]")
parserDBG.add_argument('-verbose', dest="verbosity", action="store", type=int, default=0, help="Verbosity level for DBG assembly [default: 0]")
parserDBG.add_argument("--multiple", action="store_true", help="To return the assembled sequences even if multiple solutions are found (by default, if MTG-Link returns multiple solutions, we consider 'No Assembly' as it is not possible to know which one is the correct one, except when the '--force' argument is provided)")

# IRO module options.
parserIRO = subparsers.add_parser("IRO", help="Local assembly using an Iterative Read Overlap (IRO) algorithm")
parserIRO.add_argument('-gfa', dest="gfa", action="store", help="Input GFA file (GFA 2.0) (format: xxx.gfa)", required=True)
parserIRO.add_argument('-bam', dest="bam", action="store", help="BAM file: linked-reads mapped on reference genome (format: xxx.bam)", required=True)
parserIRO.add_argument('-fastq', dest="reads", action="store", help="File (gzipped) of indexed reads (format: xxx.fastq.gz | xxx.fq.gz)", required=True)
parserIRO.add_argument('-index', dest="index", action="store", help="Barcodes index file (format: xxx.bci)", required=True)
parserIRO.add_argument('-out', dest="outDir", action="store", default="./MTG-Link_results", help="Output directory [default: ./MTG-Link_results]")
parserIRO.add_argument('-line', dest="line", action="store", type=int, help="Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]")
parserIRO.add_argument('-bxuDir', dest="bxuDir", action="store", help="Directory where the FASTQ files containing the subsample of reads are located (1 file per target) (format of FASTQ files: xxx.bxu.fastq) [to provide if the read subsampling step has already been done for this dataset]")
parserIRO.add_argument('-t', dest="threads", action="store", type=int, default=1, help="Number of threads to use for the Read Subsampling step [default: 1]")
parserIRO.add_argument('-flank', dest="flankSize", action="store", type=int, default=10000, help="Flanking sequences' size (bp) [default: 10000]")
parserIRO.add_argument('-occ', dest="minBarcOcc", action="store", type=int, default=2, help="Minimum number of occurrences in target flanking regions for a barcode to be retained in the union set [default: 2]")
parserIRO.add_argument('-ext', dest="extSize", action="store", type=int, default=500, help="Size of the extension of the target on both sides (bp); determine start/end of local assembly [default: 500]")
parserIRO.add_argument('-l', dest="maxLength", action="store", type=int, default=10000, help="Maximum assembly length (bp) (it could be a bit bigger than the length of the target to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]")
parserIRO.add_argument('-s', dest="seedSize", action="store", type=int, default=10, help="Seed size used for indexing the reads (bp) [default: 10]")
parserIRO.add_argument('-o', dest="minOverlap", action="store", type=int, default=20, help="Minimum overlapping size (bp) [default: 20]")
parserIRO.add_argument('-a', dest="abundanceMin", action="store", type=int, default=[3, 2], nargs='+', help="Minimal abundance(s) of reads used for local assembly ; extension's groups having less than this number of reads are discarded from the graph [default: [3, 2]]")
parserIRO.add_argument('-dmax', dest="maxScore", action="store", type=int, default=2, help="Maximum number of gaps/substitutions allowed in the inexact overlap between reads [default: 2]")

args = parser.parse_args()

if re.match('^.*.gfa$', args.gfa) is None:
    parser.error("\nWarning: The suffix of the GFA file should be: '.gfa'")

if re.match('^.*.bam$', args.bam) is None:
    parser.error("\nWarning: The suffix of the BAM file should be: '.bam'")


#----------------------------------------------------
# Input files and arguments
#----------------------------------------------------
# GFA 2.0 file.
gfaFile = os.path.abspath(args.gfa)
if not os.path.exists(gfaFile):
    parser.error("\nWarning: The path of the GFA file doesn't exist")
gfa_name = gfaFile.split('/')[-1]
print("\nInput GFA file: " + gfaFile)

# BAM file: linked reads mapped on current genome assembly.
bamFile = os.path.abspath(args.bam)
if not os.path.exists(bamFile): 
    parser.error("\nWarning: The path of the BAM file doesn't exist")
print("BAM file: " + bamFile)

# Reads file: file of indexed reads.
readsFile = os.path.abspath(args.reads)
if not os.path.exists(readsFile):
    parser.error("\nWarning: The path of the file of indexed reads doesn't exist")
print("File of indexed reads: " + readsFile)

# Barcodes index file.
indexFile = os.path.abspath(args.index)
if not os.path.exists(indexFile):
    parser.error("\nWarning: The path of the barcodes index file doesn't exist")
print("Barcodes index file: " + indexFile)

# Directory where the FASTQ files containing the subsample of reads are located (1 file per target).
if args.bxuDir is not None:
    bxuDir = os.path.abspath(args.bxuDir)
    if not os.path.exists(bxuDir):
        parser.error("\nWarning: The path of the directory where the FASTQ files containing the subsample of reads are located doesn't exist")
    print("Directory containing the subsample of reads: " + bxuDir)
else:
    bxuDir = ""

# Directory containing the reference sequences if any.
# if args.refDir is not None:
#     refDir = os.path.abspath(args.refDir)
#     if not os.path.exists(refDir):
#         parser.error("\nWarning: The path of the directory containing the reference sequences doesn't exist")
#     print("Directory containing the reference sequences: " + refDir)
# else:
#     refDir = ""
refDir = ""

# Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input).
if args.line is not None:
    line_gfa = args.line
else:
    line_gfa = ""


#----------------------------------------------------
# Directories for saving results
#----------------------------------------------------
cwd = os.getcwd() 

# Create the directory 'outDir', where the results will be saved.
if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)
try:
    os.chdir(args.outDir)
except OSError:
    print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()
print("\nThe results are saved in " + outDir + "\n")

# Create the subdirectory 'subsamplingDir', where the results of the read subsampling step will be saved.
subsamplingDir = outDir + "/read_subsampling"
os.mkdir(subsamplingDir)

# Create the subdirectory 'assemblyDir', where the results of the local assembly step will be saved.
assemblyDir = outDir + "/local_assembly"
os.mkdir(assemblyDir)

# Create the subdirectory 'contigDir', where the sequences of the flanking contigs will be saved.
contigDir = outDir + "/contigs"
os.mkdir(contigDir)

# Create the subdirectory 'evalDir', where the results of the qualitative evaluation step will be saved.
evalDir = outDir + "/qual_evaluation"
os.mkdir(evalDir)


#----------------------------------------------------
# Parameters
#----------------------------------------------------
try:
    module = args.MTGLink_module
    chunkSize = args.flankSize
    barcodesMinOcc = args.minBarcOcc
    extSize = args.extSize
    maxLength = args.maxLength
    threads = args.threads

    # Module DBG.
    if module == "DBG":
        minLength = args.minLength
        kmerSizeList = args.kmerSize
        abundanceThresholdList = args.abundanceThreshold
        maxNodes = args.maxNodes
        nbCores = args.nbCores
        maxMemory = args.maxMemory
        verbosity = args.verbosity

    # Module IRO.
    if module == "IRO":
        seedSize = args.seedSize
        minOverlapSize = args.minOverlap
        abundanceMinList = args.abundanceMin
        dmax = args.maxScore

except Exception as e:
    print("\nFile 'main.py': Something wrong with the parameters")
    print("Exception-")
    print(e)
    sys.exit(1)

