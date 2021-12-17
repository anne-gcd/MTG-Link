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

"""Module 'gapFilling.py': Gap Filling

The module 'gapFilling.py' enables to perform two out of the three main steps of the MTG-Link gap-filling pipeline: Local Assembly and Qualitative Evaluation.
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
    To perform the gap-filling on a specific gap. 
    This consists of two main steps: Local Assembly and Qualitative Evaluation. 
    NB: The Read Subsampling step was performed before. 

    Arg:
        - current_gap: str
            current gap identification
    
    Return/Outputs:
        - outputGFAList: list of lists
            list of lists, each list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the gap-filled sequence
    """
    #----------------------------------------------------
    # Local Assembly
    #----------------------------------------------------  
    try:
        # Perform the local assembly with the DBG (De Bruijn Graph) algorithm or the IRO (Iterative Read Overlap) algorithm. 
        ## DBG algorithm.
        if main.module == "DBG":
            gapfillingFile = localAssemblyWithDBGAlgorithm(current_gap, main.gfaFile, main.chunkSize, main.extSize, main.maxLength, main.kmerSizeList, main.abundanceThresholdList, main.maxNodes, main.nbCores, main.maxMemory, main.verbosity)

        ## IRO algorithm.
        if main.module == "IRO":
            gapfillingFile = localAssemblyWithIROAlgorithm(current_gap, main.gfaFile, main.chunkSize, main.extSize, main.maxLength, main.seedSize, main.minOverlapSize, main.abundanceMinList, main.dmax)

    except Exception as e:
        print("File 'gapFilling.py': Something wrong with the 'Local Assembly' step of the function 'fillGapByLocalAssemblyAndQualitativeEvaluation()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Qualitative Evaluation
    #----------------------------------------------------
    try:
        # If at least one solution is found, perform the Qualitative Evaluation step on the gap-filled sequence(s).
        if os.path.getsize(gapfillingFile) > 0:
            outputGFAList = qualitativeEvaluationOfTheAssembly(current_gap, main.gfaFile, main.extSize, gapfillingFile, main.module)
        
        # If no solution found, set 'outputGFAList' to empty list.
        else:
            outputGFAList = []

        # Save the current G line into the variable 'outputGFAList', only if this variable is empty.
        if len(outputGFAList) == 0:
            outputGFAList.append([str(current_gap)])

    except Exception as e:
        print("\nFile 'gapFilling.py': Something wrong with the 'Qualitative Evaluation' step of the function 'fillGapByLocalAssemblyAndQualitativeEvaluation()'")
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

