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

"""Module 'gapFilling.py': Local Assembly

The module 'gapFilling.py' enables to perform two out of the three main steps of the MTG-Link local assembly pipeline: Local Assembly and Qualitative Evaluation.
(The Read Subsampling step was performed before). 
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
import gfapy
import main
from localAssemblyDBG import localAssemblyWithDBGAlgorithm
from localAssemblyIRO import localAssemblyWithIROAlgorithm
from qualitativeEvaluation import qualitativeEvaluationOfTheAssembly


#----------------------------------------------------
# fillGapByLocalAssemblyAndQualitativeEvaluation function
#----------------------------------------------------
def fillGapByLocalAssemblyAndQualitativeEvaluation(current_gap):
    """
    To perform the local assembly on a specific gap/target. 
    This consists of two main steps: Local Assembly and Qualitative Evaluation. 
    NB: The Read Subsampling step was performed before. 

    Arg:
        - current_gap: str
            current gap/target identification
    
    Return/Outputs:
        - outputGFAList: list of lists
            list of lists, each list containing the assembled sequence's name, as well as its length, its sequence, the full solution name with the orientation sign, the k-mer size used for the DBG assembly, the beginning and ending positions of the overlap and the quality of the assembled sequence
    """
    #----------------------------------------------------
    # Local Assembly
    #----------------------------------------------------  
    try:
        # Create a list of gapfilling files 'gapfillingFilesList'.
        gapfillingFilesList = []

        # Perform the local assembly with the DBG (De Bruijn Graph) algorithm or the IRO (Iterative Read Overlap) algorithm. 
        ## DBG algorithm.
        if main.module == "DBG":
            gapfillingFile = localAssemblyWithDBGAlgorithm(current_gap, main.gfaFile, main.chunkSize, main.extSize, main.maxLength, main.minLength, main.kmerSizeList, main.abundanceThresholdList, main.maxNodes, main.nbCores, main.maxMemory, main.verbosity)
            gapfillingFilesList.append(gapfillingFile)

            # If param '--force' set by user, force search on all k-mer values provided. 
            if main.args.force:

                # Get the list of all k-mer values not already tested.
                prev_kValue = int(str(gapfillingFile).split('.bxu..insertions_filtered.fasta')[0].split('.a')[-2].split('.k')[-1])
                updated_kmerSizeList = [k for k in main.kmerSizeList if k < prev_kValue]
            
                # Continue to search for an assembly sequence for all k-mer values not already tested.
                while len(updated_kmerSizeList) != 0:

                    # Local assembly for the remaining k-mer values.
                    gapfillingFile = localAssemblyWithDBGAlgorithm(current_gap, main.gfaFile, main.chunkSize, main.extSize, main.maxLength, main.minLength, updated_kmerSizeList, main.abundanceThresholdList, main.maxNodes, main.nbCores, main.maxMemory, main.verbosity)
                    gapfillingFilesList.append(gapfillingFile)

                    # Get the list of all k-mer values not already tested.
                    prev_kValue = int(str(gapfillingFile).split('.bxu..insertions_filtered.fasta')[0].split('.a')[-2].split('.k')[-1])
                    updated_kmerSizeList = [k for k in main.kmerSizeList if k < prev_kValue]

        ## IRO algorithm.
        if main.module == "IRO":
            gapfillingFile = localAssemblyWithIROAlgorithm(current_gap, main.gfaFile, main.chunkSize, main.extSize, main.maxLength, main.seedSize, main.minOverlapSize, main.abundanceMinList, main.dmax)
            gapfillingFilesList.append(gapfillingFile)

    except Exception as e:
        print("File 'gapFilling.py': Something wrong with the 'Local Assembly' step of the function 'fillGapByLocalAssemblyAndQualitativeEvaluation()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Qualitative Evaluation
    #----------------------------------------------------
    try:
        # Create an empty list 'outputGFAList'.
        ## List of lists, each list containing the assembled sequence's name, as well as its length, its sequence, the full solution name with the orientation sign, the k-mer size used for the DBG assembly, the beginning and ending positions of the overlap and the quality of the assembled sequence.
        outputGFAList = []

        # Iterate over the 'gapfillingFile' files.
        for gapfillingFile in gapfillingFilesList:

            # If at least one solution is found, perform the Qualitative Evaluation step on the assembled sequence(s).
            if os.path.getsize(gapfillingFile) > 0:
                outputGFAList = qualitativeEvaluationOfTheAssembly(current_gap, main.gfaFile, main.extSize, gapfillingFile, main.module, outputGFAList)

        # Save the current G line into the variable 'outputGFAList', only if this variable is empty.
        if len(outputGFAList) == 0:
            outputGFAList.append([str(current_gap)])

    except Exception as e:
        print("File 'gapFilling.py': Something wrong with the 'Qualitative Evaluation' step of the function 'fillGapByLocalAssemblyAndQualitativeEvaluation()'")
        print("Exception-")
        print(e)
        sys.exit(1)
        
    # Go in the 'outDir' directory.
    try:
        os.chdir(main.outDir)
    except OSError as err:
            print("File 'gapFilling.py': Something wrong with specified directory 'outDir'. \nOSError-{}".format(err))
            sys.exit(1)


    return outputGFAList

