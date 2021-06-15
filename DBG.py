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

"""Module 'DBG.py': De Bruijn Graph (DBG) algorithm

The module 'DBG.py' enables to perform the local assembly step of the MTG-Link gap-filling pipeline, using a De Bruijn Graph (DBG) algorithm. 
The DBG algorithm is performed with the 'fill' module of the software MindTheGap, using the subsample of linked reads obtained during the first step of the MTG-Link pipeline.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
from gfapy.sequence import rc
from Bio import SeqIO
from main import gfa_name, subsamplingDir, assemblyDir, chunk_size, ext_size, kmer_sizeList, abundance_thresholdList, max_nodes, nb_cores, max_memory, verbosity


#----------------------------------------------------
# mtg_fill function
#----------------------------------------------------
def mtg_fill(gap_label, reads_file, bkpt_file, k, a, max_nodes, max_length, nb_cores, max_memory, verbosity, output_prefix):
    """
    To execute the MindTheGap fill module, that relies on a De Bruijn graph data structure to represent the input read sequences.
    The local assembly step is performed between the sequences surrounding the extended gap (e.g. the kmers of the breakpoint file), in both orientations.

    Args:
        - gap_label: str
            label of the gap
        - reads_file: file
            reads of the union
        - bkpt_file: file
            breakpoint sequences, with offset of size k removed
        - k: int
            k-mer size
        - a: int
            minimal abundance threshold for solide k-mers
        - max_nodes: int
            maximum number of nodes in contig graph
        - max_length: int
            maximum assembly length (bp)
        - nb_cores: int
            number of cores
        - max_memory: int
            maximum memory for graph building (in MBytes)
        - verbosity: int
            verbosity level
        - output_prefix: str
            this function will save the assembly results on several files whose prefix is the 'output_prefix' name

    Output:
        - output_prefix: str
            prefix of files containing the assembly results obtained with `MindTheGap fill`

    Return:
        - res: str
            - the file containing the assembled sequence(s) and a Boolean value equal to True
              if a solution is found (e.g. we arrived to stop kmer)
            OR
            - the reason why the gap-filling failed and a Boolean value equal to False
              if no solution is found
    """
    try:
        # `MindTheGap fill`.
        ## The option '-fwd-only' is used to avoid redundancies as MTG-Link already perform the local assembly step in both orientations
        if max_memory == 0:
            command = ["/home/genouest/inra_umr1349/aguichard/Gapfilling/DBG/scripts/MindTheGap/build/bin/MindTheGap", "fill", "-in", reads_file, "-bkpt", bkpt_file, \
                        "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(max_nodes), "-max-length", str(max_length), \
                        "-nb-cores", str(nb_cores), "-verbose", str(verbosity), "-fwd-only", "-out", output_prefix]
        else:
            command = ["/home/genouest/inra_umr1349/aguichard/Gapfilling/DBG/scripts/MindTheGap/build/bin/MindTheGap", "fill", "-in", reads_file, "-bkpt", bkpt_file, \
                        "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(max_nodes), "-max-length", str(max_length), \
                        "-nb-cores", str(nb_cores), "-max-memory", str(max_memory), "-verbose", str(verbosity), "-fwd-only", "-out", output_prefix]
        mtgfillLog = str(gap_label) + "_mtgfill.log"

        with open(mtgfillLog, "a") as log:
            subprocess.run(command, stderr=log)
            output = subprocess.check_output(command)

        # Remove the raw files obtained from `MindTheGap fill`.
        if os.path.getsize(mtgfillLog) <= 0:
            subprocess.run(["rm", mtgfillLog])

        # If we find a complete gap-filled sequence (e.g. we reach the kmer stop), return the assembly sequence along with True.
        if os.path.getsize(assemblyDir +"/"+ output + ".insertions.fasta") > 0:
            res = os.path.abspath(assemblyDir +"/"+ output + ".insertions.fasta")
            return res, True

        # If we don't find a complete gap-filled sequence, return "Gap-filling not completed..." along with False.
        else:
            res = "Gap-filling not completed..."
            return res, False

    except Exception as e:
        print("\nFile 'DBG.py': Something wrong with the function 'mtg_fill()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# dbg_assembly function
#----------------------------------------------------
def dbg_assembly(gap_label, gap, left_scaffold, right_scaffold, seq_L, seq_R, max_length):
    """
    To perform the Local Assembly step using a De Bruijn Graph (DBG) algorithm.  
    The DBG algorithm is performed with the 'fill' module of the software MindTheGap. This module is executed on the reads of the union, in breakpoint mode.
    This consists of two main steps: get the Breakpoing file, perform the local assembly with `MindTheGap fill`.

    Args:
        - gap_label: str
            label of the gap
        - gap: object
            Gap object for the current gap
        - left_scaffold: object
            Scaffold object for the left scaffold
        - right_scaffold: object
            Scaffold object for the right scaffold
        - seq_L: str
            left flanking contig sequence
        - seq_R: str
            right flanking contig sequence
        - max_length: int
            maximum assembly length (bp)

    Return:
        - gapfillingFile: file
            file containing the obtained gap-filled sequence(s)
    """
    # Iterate over the kmer values, starting with the highest.
    for k in kmer_sizeList:
    
        #----------------------------------------------------
        # Breakpoint file, with offset of size k removed
        #----------------------------------------------------
        try:
            # Get a breakpoint file containing the input sequences (start and stop kmers).
            bkptFile = "{}.{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, str(gap_label), gap.length, chunk_size, k)
            with open(bkptFile, "w") as bkpt:

                # Left kmer and Reverse Right kmer (dependent on the orientation of the left scaffold).
                line1 = ">bkpt1_GapID.{}_Gaplen.{} left_kmer.{}_len.{} offset_rm\n".format(str(gap_label), gap.length, left_scaffold.name, k)
                line2 = seq_L[(left_scaffold.slen - ext_size - k):(left_scaffold.slen - ext_size)]
                line7 = "\n>bkpt2_GapID.{}_Gaplen.{} right_kmer.{}_len.{} offset_rm\n".format(str(gap_label), gap.length, left_scaffold.name, k)
                line8 = str(rc(seq_L)[ext_size:(ext_size + k)])

                # Right kmer and Reverse Left kmer (dependent on the orientation of the right scaffold).
                line3 = "\n>bkpt1_GapID.{}_Gaplen.{} right_kmer.{}_len.{} offset_rm\n".format(str(gap_label), gap.length, right_scaffold.name, k)
                line4 = seq_R[ext_size:(ext_size + k)]
                line5 = "\n>bkpt2_GapID.{}_Gaplen.{} left_kmer.{}_len.{} offset_rm\n".format(str(gap_label), gap.length, right_scaffold.name, k)
                line6 = str(rc(seq_R)[(right_scaffold.slen - ext_size - k):(right_scaffold.slen - ext_size)])

                bkpt.writelines([line1, line2, line3, line4, line5, line6, line7, line8])

        except Exception as e:
            print("\nFile 'DBG.py': Something wrong with the input bkpt file (containing start/stop k-mers) creation")
            print("Exception-")
            print(e)
            sys.exit(1)


        #----------------------------------------------------
        # `MindTheGap fill`
        #----------------------------------------------------
        # Iterate over the abundance threshold values, starting with the highest.
        for a in abundance_thresholdList:

            try:
                print("\nGapfilling of {} for k={} and a={} (union)".format(str(gap_label), k, a))
                
                # Input and output files.
                union_readsFile = "{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, str(gap_label), gap.length, chunk_size)
                subreadsFile = os.path.join(subsamplingDir, union_readsFile)
                output_prefix = "{}.{}.g{}.c{}.k{}.a{}.bxu".format(gfa_name, str(gap_label), gap.length, chunk_size, k, a)

                # Determine the maximum assembly length (bp) if the gap length is known.
                ## Add twice the extension size to the gap length (extension on both sides of the gap) and twice the length of reads (e.g. 2x 150 bp) to be large
                if gap.length >= max_length:
                    max_length = gap.length + 2*ext_size + 2*150
                
            except Exception as e:
                print("\nFile 'DBG.py': Something wrong with the initialization of the DBG algorithm (e.g. initialization of `MindTheGap fill`)")
                print("Exception-")
                print(e)
                sys.exit(1)
            
            # Perform the local assembly step with `MindTheGap fill`.
            res, success = mtg_fill(gap_label, subreadsFile, bkptFile, k, a, max_nodes, max_length, nb_cores, max_memory, verbosity, output_prefix)

            # Case of successful gap-filling.
            if success:
                break

        # Case of successful gap-filling.
        if success:
            break
                    

    #----------------------------------------------------
    # Output from `MindTheGap fill`
    #----------------------------------------------------
    try:
        # Output file containing the gap-filled sequence(s)
        insertion_file = "{}.{}.g{}.c{}.k{}.a{}.bxu.insertions.fasta".format(gfa_name, str(gap_label), gap.length, chunk_size, k, a)

        # Case of unsuccessful gap-filling.
        if not success:
            print("\n{}: {}".format(gap_label, res))
            gapfillingFile = insertion_file

        # Case of successful gap-filling.
        if success:
            print("\n{}: Successful Gap-filling !". format(gap_label))

            # Save and pre-process the file containing the gap-filled sequence(s) for further qualitative evaluation.
            ## Modify the 'insertion_file' and save it to a new file ('gapfillingFile') so that the 'solution x/y' part appears in record.id (and not just in record.description)
            gapfillingFile = os.path.abspath(assemblyDir +"/"+ output_prefix + "..insertions.fasta")
            with open(insertion_file, "r") as original, open(gapfillingFile, "w") as corrected:
                records = SeqIO.parse(original, "fasta")
                for record in records:
                    if "solution" in record.description:
                        record.id = record.id + "_sol_" + record.description.split(" ")[-1]
                    else:
                        record.id = record.id + "_sol_1/1"
                    SeqIO.write(record, corrected, "fasta")
    
    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the processing of the output from the IRO algorithm")
        print("Exception-")
        print(e)
        sys.exit(1)


    return gapfillingFile


