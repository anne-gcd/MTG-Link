#!/usr/bin/env python3
#*****************************************************************************
#  Name: MTG-Link
#  Description: Gap-filling tool for draft genome assemblies, dedicated to 
#  linked read data.
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

MTGLINK_VERSION = "1.2.0"

parser = argparse.ArgumentParser(prog="mtglink.py", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=("Gapfilling with linked read data, using either a De Bruijn Graph (DBG) or an Iterative Read Overlap (IRO) algorithm"))

parser.add_argument('-v', action="version", version='%(prog)s {version}'.format(version=MTGLINK_VERSION))

subparsers = parser.add_subparsers(dest="MTGLink_module", help="MTGLink_module used for the Local Assembly step", required=True)

# DBG module options.
parserDBG = subparsers.add_parser("DBG", help="Gap-filling using a De Bruijn Graph (DBG) algorithm")
parserDBG.add_argument('-gfa', dest="gfa", action="store", help="Input GFA file (GFA 2.0) (format: xxx.gfa)", required=True)
parserDBG.add_argument('-c', dest="chunk", action="store", type=int, default=5000, help="Chunk size (bp) [default: 5000]")
parserDBG.add_argument('-bam', dest="bam", action="store", help="BAM file: linked reads mapped on current genome assembly (format: xxx.bam)", required=True)
parserDBG.add_argument('-fastq', dest="reads", action="store", help="File of indexed reads (format: xxx.fastq | xxx.fq)", required=True)
parserDBG.add_argument('-index', dest="index", action="store", help="Barcodes index file (format: xxx.shelve)", required=True)
parserDBG.add_argument('-f', dest="freq", action="store", type=int, default=2, help="Minimal frequence of barcodes observed in the union set from the two flanking gap sequences [default: 2]")
parserDBG.add_argument('-out', dest="outDir", action="store", default="./mtglink_results", help="Output directory [default: ./mtglink_results]")
parserDBG.add_argument('-refDir', dest="refDir", action="store", help="Directory containing the reference sequences if any [optional]")
parserDBG.add_argument('-line', dest="line", action="store", type=int, help="Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]")
parserDBG.add_argument('-rbxu', dest="rbxu", action="store", help="File containing the reads of the union of the corresponding gap (if already extracted) [optional]")
parserDBG.add_argument('-ext', dest="extension", action="store", type=int, default=500, help="Size of the extension of the gap on both sides (bp); determine start/end of gapfilling [default: 500]")
parserDBG.add_argument('-l', dest="max_length", action="store", type=int, default=10000, help="Maximum assembly length (bp) (it could be a bit bigger than the length of the gap to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]")
parserDBG.add_argument('-k', dest="kmer_size", action="store", type=int, default=[51, 41, 31, 21],  nargs='+', help="k-mer size(s) used for gap-filling [default: [51, 41, 31, 21]]")
parserDBG.add_argument('-a', dest="abundance_threshold", action="store", type=int, default=[3, 2], nargs='+', help="Minimal abundance threshold for solid k-mers [default: [3, 2]]")
parserDBG.add_argument("--force", action="store_true", help="To force search on all '-k' values provided")
parserDBG.add_argument('-max-nodes', dest="max_nodes", action="store", type=int, default=1000, help="Maximum number of nodes in contig graph [default: 1000]")
parserDBG.add_argument('-nb-cores', dest="nb_cores", action="store", type=int, default=1, help="Number of cores for DBG assembly [default: 1]")
parserDBG.add_argument('-max-memory', dest="max_memory", action="store", type=int, default=0, help="Maximum memory for graph building (in MBytes) [default: 0]")
parserDBG.add_argument('-verbose', dest="verbosity", action="store", type=int, default=0, help="Verbosity level for DBG assembly [default: 0]")

# IRO module options.
parserIRO = subparsers.add_parser("IRO", help="Gap-filling using an Iterative Read Overlap (IRO) algorithm")
parserIRO.add_argument('-gfa', dest="gfa", action="store", help="Input GFA file (GFA 2.0) (format: xxx.gfa)", required=True)
parserIRO.add_argument('-c', dest="chunk", action="store", type=int, default=5000, help="Chunk size (bp) [default: 5000]")
parserIRO.add_argument('-bam', dest="bam", action="store", help="BAM file: linked reads mapped on current genome assembly (format: xxx.bam)", required=True)
parserIRO.add_argument('-fastq', dest="reads", action="store", help="File of indexed reads (format: xxx.fastq | xxx.fq)", required=True)
parserIRO.add_argument('-index', dest="index", action="store", help="Barcodes index file (format: xxx.shelve)", required=True)
parserIRO.add_argument('-f', dest="freq", action="store", type=int, default=2, help="Minimal frequence of barcodes observed in the union set from the two flanking gap sequences [default: 2]")
parserIRO.add_argument('-out', dest="outDir", action="store", default="./mtglink_results", help="Output directory [default: ./mtglink_results]")
parserIRO.add_argument('-refDir', dest="refDir", action="store", help="Directory containing the reference sequences if any [optional]")
parserIRO.add_argument('-line', dest="line", action="store", type=int, help="Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]")
parserIRO.add_argument('-rbxu', dest="rbxu", action="store", help="File containing the reads of the union of the corresponding gap (if already extracted) [optional]")
parserIRO.add_argument('-ext', dest="extension", action="store", type=int, default=500, help="Size of the extension of the gap on both sides (bp); determine start/end of gapfilling [default: 500]")
parserIRO.add_argument('-l', dest="max_length", action="store", type=int, default=10000, help="Maximum assembly length (bp) (it could be a bit bigger than the length of the gap to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]")
parserIRO.add_argument('-s', dest="seed_size", action="store", type=int, default=10, help="Seed size used for indexing the reads (bp) [default: 10]")
parserIRO.add_argument('-o', dest="min_overlap", action="store", type=int, default=20, help="Minimum overlapping size (bp) [default: 20]")
parserIRO.add_argument('-a', dest="abundance_min", action="store", type=int, default=[3, 2], nargs='+', help="Minimal abundance(s) of reads used for gapfilling ; extension's groups having less than this number of reads are discarded from the graph [default: [3, 2]]")
parserIRO.add_argument('-dmax', dest="max_score", action="store", type=int, default=2, help="Maximum number of gaps/substitutions allowed in the inexact overlap between reads [default: 2]")

args = parser.parse_args()

if re.match('^.*.gfa$', args.input_gfa) is None:
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

# Directory containing the reference sequences if any.
if args.refDir is not None:
    refDir = os.path.abspath(args.refDir)
    if not os.path.exists(refDir):
        parser.error("\nWarning: The path of the directory containing the reference sequences doesn't exist")
    print("Directory containing the reference sequences: " + refDir)
else:
    refDir = ""

# File containing the reads of the union of the corresponding gap (if already extracted).
if args.rbxu is not None:
    rbxuFile = os.path.abspath(args.rbxu)
    if not os.path.exists(rbxuFile):
        parser.error("\nWarning: The path of the file containing the reads of the union of the corresponding gap (already extracted) doesn't exist")
    print("File containing the reads of the union of the corresponding gap (already extracted): " + rbxuFile)
else:
    rbxuFile = ""

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
    chunk_size = args.chunk
    barcodes_min_freq = args.freq
    ext_size = args.extension
    max_length = args.max_length
    kmer_sizeList = args.kmer_size
    abundance_thresholdList = args.abundance_threshold
    max_nodes = args.max_nodes
    nb_cores = args.nb_cores
    max_memory = args.max_memory
    verbosity = args.verbosity
    seed_size = args.seed_size
    min_overlap = args.min_overlap
    abundance_minList = args.abundance_min
    dmax = args.max_score

except Exception as e:
    print("\nFile 'main.py': Something wrong with the parameters")
    print("Exception-")
    print(e)
    sys.exit(1)

