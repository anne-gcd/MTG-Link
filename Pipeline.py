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

"""Module 'Pipeline.py': Pipeline

The module 'Pipeline.py' enables to perform main steps of the MTG-Link gap-filling pipeline.
This consists of three main steps: Read Subsampling, Local Assembly, Qualitative Evaluation. 
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
import gfapy
from main import gfaFile, outDir, subsamplingDir, assemblyDir, contigDir, module, chunk_size, max_length
from helpers import Gap, Scaffold
from ReadSubsampling import read_subsampling
from DBG import dbg_assembly
from IRO import iro_assembly
from QualEval import qual_eval


#----------------------------------------------------
# gapfilling function
#----------------------------------------------------
def gapfilling(current_gap):
    """
    To perform the gap-filling on a specific gap. 
    This consists of three main steps: Read Subsampling, Local Assembly, Qualitative Evaluation. 

    Arg:
        - current_gap: str
            current gap identification
    
    Return/Outputs:
        - union_summary: list
            list containing the gap ID, the names of the left and right flanking sequences, the gap size, the chunk size, and the number of barcodes and reads extracted on the chunks to perform the gap-filling
        - output_for_gfa: list
            list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the sequence
    """
    #----------------------------------------------------
    # Pre-processing
    #----------------------------------------------------
    try:
        os.chdir(outDir)
    except OSError:
        print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
        sys.exit(1)
    
    try:
        # Open the input GFA file to get the corresponding Gap line ('G' line).
        gfa = gfapy.Gfa.from_file(gfaFile)
        for _gap_ in gfa.gaps:
            if str(_gap_) == current_gap:
                current_gap = _gap_
                ## Create the object 'gap' from the class 'Gap'
                gap = Gap(current_gap)

        # Get some information on the current gap we are working on.
        gap.info()
        gap_label = gap.label()

        # Create two objects ('left_scaffold' and 'right_scaffold') from the class 'Scaffold'.
        left_scaffold = Scaffold(current_gap, gap.left, gfaFile)
        right_scaffold = Scaffold(current_gap, gap.right, gfaFile)

        # If chunk size larger than length of scaffold(s), set the chunk size to the minimal scaffold length.
        ## Left chunk
        if chunk_size > left_scaffold.slen:
            print("\nWarning for {}: The chunk size you provided is higher than the length of the left scaffold. Thus, for the left scaffold, the barcodes will be extracted on its whole length".format(gap_label))
            chunk_L = left_scaffold.slen
        else:
            chunk_L = chunk_size
        ## Right chunk
        if chunk_size > right_scaffold.slen:
            print("\nWarning for {}: The chunk size you provided is higher than the length of the right scaffold. Thus, for the right scaffold, the barcodes will be extracted on its whole length".format(gap_label))
            chunk_R = right_scaffold.slen
        else:
            chunk_R = chunk_size

    except Exception as e:
        print("\nFile 'Pipeline.py': Something wrong with the Pre-processing step")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Read Subsampling
    #----------------------------------------------------
    try:
        os.chdir(subsamplingDir)
    except OSError:
        print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
        sys.exit(1)

    try:
        # Get the list 'union_summary' containing a summary of the Read Subsampling step.
        union_summary = read_subsampling(gap_label, gap, left_scaffold, right_scaffold, chunk_L, chunk_R)
    
    except Exception as e:
        print("\nFile 'Pipeline.py': Something wrong with the Read Subsampling step")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Local Assembly
    #----------------------------------------------------  
    try:
        os.chdir(assemblyDir)
    except OSError:
        print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
        sys.exit(1)      

    try:
        # Get the flanking contigs sequences.
        seq_L = str(left_scaffold.sequence())
        seq_R = str(right_scaffold.sequence())

        # Determine which module to execute: DBG or IRO.
        if module == "DBG":
            gapfillingFile = dbg_assembly(gap_label, gap, left_scaffold, right_scaffold, seq_L, seq_R, max_length)

        if module == "IRO":
            gapfillingFile = iro_assembly(gap_label, gap, left_scaffold, right_scaffold, seq_L, seq_R, max_length)

    except Exception as e:
        print("\nFile 'Pipeline.py': Something wrong with the Local Assembly step")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Qualitative Evaluation
    #----------------------------------------------------
    try:
        # If at least one solution is found, perform the Qualitative Evaluation step on the gap-filled sequence(s).
        if os.path.getsize(gapfillingFile) > 0:
            output_for_gfa = qual_eval(gap_label, gap, left_scaffold, right_scaffold, seq_L, seq_R, gapfillingFile)
        
        # If no solution found, set 'output_for_gfa' to empty list.
        else:
            output_for_gfa = []

        # Save the current G line into the variable 'output_for_gfa' only if this variable is empty.
        if len(output_for_gfa) == 0:
            output_for_gfa.append([str(current_gap)])

        # Change directory.
        try:
            os.chdir(outDir)
        except OSError:
            print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
            sys.exit(1)

    except Exception as e:
        print("\nFile 'Pipeline.py': Something wrong with the Qualitative Evaluation step")
        print("Exception-")
        print(e)
        sys.exit(1)
        
    
    return union_summary, output_for_gfa

