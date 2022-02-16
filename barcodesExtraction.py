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

"""Module 'barcodesExtraction.py': Barcodes Extraction

The module 'barcodesExtraction.py' enables to extract the barcodes observed in chunk regions surrounding the gap.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
import gfapy
from helpers import Gap, Scaffold


#----------------------------------------------------
# extractBarcodesWithLRezExtract function
#----------------------------------------------------
def extractBarcodesWithLRezExtract(bam, gapLabel, region, barcodesOccurrencesDict):
    """
    To extract the barcodes of reads mapping on chunk regions, with `LRez extract`. 
    `LRez extract` enables to extract the list of barcodes from a given region of a BAM file.

    Args:
        - bam: file
            indexed BAM file obtained after mapping the linked reads onto the draft assembly
        - gapLabel: str
            label of the gap
        - region: str
            chunk region on which to extract the barcodes
        - barcodesOccurrencesDict: dict
            this function will output a dictionary 'barcodesOccurrencesDict' containing the occurences of each barcode
            key = barcode sequence ; value = number of occurences

    Return/Output:
        - barcodesOccurrencesDict: dict
            dictionary 'barcodesOccurrencesDict' containing the occurences for each barcode extracted on the chunk region
            key = barcode sequence ; value = number of occurences
    """
    try:
        # LRez extract. 
        ## Parameter '-d' to include duplicate barcodes.
        command = ["LRez", "extract", "--bam", bam, "--region", region, "-d"]
        extractLog = str(gapLabel) + "_LRezExtract.log"
        tmpBarcodesFile = str(gapLabel) + "_LRezExtract_stdout.txt"

        try:
            with open(tmpBarcodesFile, "w+") as f, open(extractLog, "a") as log:
                subprocess.run(command, stdout=f, stderr=log)
                f.seek(0)

                # Save the barcodes and their occurences in the dict 'barcodesOccurrencesDict'.
                for line in f.readlines():
                    ## Remove the '\n' at the end of the sequence (if any).
                    barcodeSeq = line.split('\n')[0]
                    ## Count occurrences of each barcode and save them in the dict 'barcodesOccurrencesDict'.
                    if barcodeSeq in barcodesOccurrencesDict:
                        barcodesOccurrencesDict[barcodeSeq] += 1
                    else:
                        barcodesOccurrencesDict[barcodeSeq] = 1

        except IOError as err:
            print("File 'barcodesExtraction.py', function 'extractBarcodesWithLRezExtract()': Unable to open or write to the files {} and/or {}. \nIOError-{}".format(str(tmpBarcodesFile), str(extractLog), err))
            sys.exit(1)

        # Remove the raw files obtained from `LRez extract`. 
        subprocess.run(["rm", tmpBarcodesFile])
        if os.path.getsize(extractLog) == 0:
            subprocess.run(["rm", extractLog])

        return barcodesOccurrencesDict

    except Exception as e:
        print("File 'barcodesExtraction.py': Something wrong with the function 'extractBarcodesWithLRezExtract()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# extractBarcodesFromChunkRegions function
#----------------------------------------------------
def extractBarcodesFromChunkRegions(current_gap, gfaFile, bamFile, chunkSize, barcodesMinFreq):
    """
    To extract the barcodes of reads mapping on chunk regions. 

    Args:
        - current_gap: str
            current gap identification
        - gfaFile: file
            GFA file containing the gaps' coordinates
        - bamFile: file
            indexed BAM file obtained after mapping the linked reads onto the draft assembly
        - chunkSize: int
            size of the chunk region
        - barcodesMinFreq: int
            minimal frequence of barcodes observed in the union set from the two flanking gap sequences

    Return:
        - unionBarcodesFile: file
            file containing the extracted barcodes of the union of both left and right gap flanking sequences
    """
    #----------------------------------------------------
    # Pre-Processing
    #----------------------------------------------------
    try:
        # Open the input GFA file.
        gfa = gfapy.Gfa.from_file(gfaFile)
        if not gfa:
            print("File 'barcodesExtraction.py', function 'extractBarcodesFromChunkRegions()': Unable to open the input GFA file {}.".format(str(gfaFile)), file=sys.stderr)
            sys.exit(1)

        # Get the corresponding Gap line ('G' line).
        for _gap_ in gfa.gaps:
            if str(_gap_) == current_gap:
                current_gap = _gap_
                ## Create the object 'gap' from the class 'Gap'
                gap = Gap(current_gap)
                if not gap:
                    print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Unable to create the object 'gap' from the class 'Gap'.", file=sys.stderr)
                    sys.exit(1)

        # Get some information on the current gap we are working on.
        gap.info()
        gapLabel = gap.label()

        # Create two objects ('leftScaffold' and 'rightScaffold') from the class 'Scaffold'.
        leftScaffold = Scaffold(current_gap, gap.left, gfaFile)
        if not leftScaffold:
            print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Unable to create the object 'leftScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)
        rightScaffold = Scaffold(current_gap, gap.right, gfaFile)
        if not rightScaffold:
            print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Unable to create the object 'rightScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)

        # If chunk size larger than length of scaffold(s), set the chunk size to the minimal scaffold length.
        ## Left chunk
        if chunkSize > leftScaffold.slen:
            print("Warning for {}: The chunk size you provided is higher than the length of the left scaffold. Thus, for the left scaffold, the barcodes will be extracted on its whole length".format(gapLabel))
            chunk_L = leftScaffold.slen
        else:
            chunk_L = chunkSize
        ## Right chunk
        if chunkSize > rightScaffold.slen:
            print("Warning for {}: The chunk size you provided is higher than the length of the right scaffold. Thus, for the right scaffold, the barcodes will be extracted on its whole length".format(gapLabel))
            chunk_R = rightScaffold.slen
        else:
            chunk_R = chunkSize

    except Exception as e:
        print("File 'barcodesExtraction.py': Something wrong with the 'Pre-Processing' step of the function 'extractBarcodesFromChunkRegions()'")
        print("Exception-{}".format(e))
        sys.exit(1)

    
    #----------------------------------------------------
    # Extract Barcodes
    #----------------------------------------------------    
    try:
        # Initiate a dictionary to count the occurences of each barcode extracted on the chunk regions.
        barcodesOccurrencesDict = {}
        
        # Obtain the left barcodes extracted on the left region (left chunk) and store the barcodes and their occurences in the dict 'barcodes_occ'.
        leftRegion = leftScaffold.chunk(chunk_L)
        if not leftRegion:
            print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Unable to obtain the left region (left chunk).", file=sys.stderr)
            sys.exit(1)
        extractBarcodesWithLRezExtract(bamFile, gapLabel, leftRegion, barcodesOccurrencesDict)

        # Obtain the right barcodes extracted on the right region (right chunk) and store the barcodes and their occurences in the dict 'barcodes_occ'.
        rightRegion = rightScaffold.chunk(chunk_R)
        if not rightRegion:
            print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Unable to obtain the right region (right chunk).", file=sys.stderr)
            sys.exit(1)
        extractBarcodesWithLRezExtract(bamFile, gapLabel, rightRegion, barcodesOccurrencesDict)

        if len(barcodesOccurrencesDict) == 0:
            print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Error while extracting the barcodes with `LRez extract`.", file=sys.stderr)
            sys.exit(1)

        # Do the union of the barcodes on both left and right regions (e.g. both left and right chunks). 
        gfa_name = gfaFile.split('/')[-1]
        unionBarcodesFile = "{}.{}.g{}.c{}.f{}.bxu".format(gfa_name, str(gapLabel), gap.length, chunkSize, barcodesMinFreq)
        try:
            with open(unionBarcodesFile, "w") as unionBarcFile:
                ## Filter barcodes by the minimal frequence of barcodes observed in the union set from the two flanking gap sequences ('barcodesMinFreq')
                for (barcode, occurrence) in barcodesOccurrencesDict.items():
                    if occurrence >= barcodesMinFreq:
                        unionBarcFile.write(barcode + "\n")
        except IOError as err:
            print("File 'barcodesExtraction.py, function 'extractBarcodesFromChunkRegions()': Unable to open or write to the output 'unionBarcodesFile' {}. \nIOError-{}".format(str(unionBarcodesFile), err))
            sys.exit(1)

    except Exception as e:
        print("File 'barcodesExtraction.py': Something wrong with the 'Extract Barcodes' step of the function 'extractBarcodesFromChunkRegions()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    return unionBarcodesFile

