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

"""Module 'localAssemblyDBG.py': Local Assembly With the DBG (De Bruijn Graph) algorithm

The module 'localAssemblyDBG.py' enables to perform the local assembly step of the MTG-Link local assembly pipeline, using a DBG (De Bruijn Graph) algorithm. 
The DBG algorithm is performed with the 'fill' module of the software MindTheGap, using the subsample of linked reads obtained during the 'Read Subsampling' step of the MTG-Link pipeline.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
from turtle import right
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO
import main
from helpers import Gap, Scaffold, getMostRepresentedKmer


#----------------------------------------------------
# fillGapWithMTGFill function
#----------------------------------------------------
def fillGapWithMTGFill(gapLabel, readsFile, bkptFile, k, a, maxLength, maxNodes, nbCores, maxMemory, verbosity, outputPrefix):
    """
    To execute the MindTheGap fill module, that relies on a De Bruijn graph data structure to represent the input read sequences.
    The local assembly step is performed between the sequences surrounding the extended gap/target (e.g. the kmers of the breakpoint file), in both orientations.

    Args:
        - gapLabel: str
            label of the gap/target
        - readsFile: file
            reads file used for the local assembly step with `MindTheGap fill`
        - bkptFile: file
            breakpoint file containing the input sequences for the local assembly with `MindTheGap fill` (start and stop kmers)
            (breakpoint sequences with offset of size k removed)
        - k: int
            k-mer size
        - a: int
            minimal abundance threshold for solid k-mers
        - maxLength: int
            maximum assembly length (bp)
        - maxNodes: int
            maximum number of nodes in contig graph
        - nbCores: int
            number of cores
        - maxMemory: int
            maximum memory for graph building (in MBytes)
        - verbosity: int
            verbosity level
        - outputPrefix: str
            prefix of the output files obtained from `MindTheGap fill`. 

    Return:
        - res: str
            - the file containing the assembled sequence(s) (and so the Boolean value equal to True)
              if a solution is found (e.g. we arrived to stop kmer)
            OR
            - the sentence "Local Assembly not completed..." (and so the Boolean value equal to False)
              if no solution is found
        - True/False: boolean
    """
    try:
        # MindTheGap fill.
        ## The option '-fwd-only' is used to avoid redundancies as MTG-Link already performs the local assembly step in both orientations.
        if maxMemory == 0:
            command = ["MindTheGap", "fill", "-in", str(readsFile), "-bkpt", str(bkptFile), \
                        "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(maxNodes), "-max-length", str(maxLength), \
                        "-nb-cores", str(nbCores), "-verbose", str(verbosity), "-fwd-only", "-out", outputPrefix]
        else:
            command = ["MindTheGap", "fill", "-in", str(readsFile), "-bkpt", str(bkptFile), \
                        "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(maxNodes), "-max-length", str(maxLength), \
                        "-nb-cores", str(nbCores), "-max-memory", str(maxMemory), "-verbose", str(verbosity), "-fwd-only", "-out", outputPrefix]
        mtgFillLog = str(gapLabel) + "_MTGFill.log"

        try:
            with open(mtgFillLog, "a") as log:
                subprocess.run(command, stderr=log)
                #output = subprocess.check_output(command)
        except IOError as err:
            print("File 'localAssemblyDBG.py', function 'fillGapWithMTGFill()': Unable to open or write to the file {}. \nIOError-{}".format(str(mtgFillLog), err))
            sys.exit(1)

        # Remove the raw files obtained from `MindTheGap fill`.
        if os.path.getsize(mtgFillLog) <= 0:
            subprocess.run(["rm", mtgFillLog])

        # If we find a complete assembled sequence (e.g. we reach the kmer stop), return the assembly sequence along with True.
        if (os.path.getsize(main.assemblyDir +"/"+ outputPrefix + ".insertions.fasta") > 0):
            res = os.path.abspath(main.assemblyDir +"/"+ outputPrefix + ".insertions.fasta")
            return res, True

        # If we don't find a complete assembled sequence, return "Local Assembly not completed..." along with False.
        else:
            res = "Local Assembly not completed..."
            return res, False

    except Exception as e:
        print("File 'localAssemblyDBG.py': Something wrong with the function 'fillGapWithMTGFill()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# localAssemblyWithDBGAlgorithm function
#----------------------------------------------------
def localAssemblyWithDBGAlgorithm(current_gap, gfaFile, chunkSize, extSize, maxLength, minLength, kmerSizeList, abundanceThresholdList, maxNodes, nbCores, maxMemory, verbosity):
    """
    To perform the Local Assembly step using a DBG (De Bruijn Graph) algorithm.
    The DBG algorithm is performed with the 'fill' module of the software MindTheGap. This module is executed on the subsample of reads retrieved during the 'Read Subsampling' step, in breakpoint mode.
    This consists of four main steps: Pre-processing of the current gap/target, getting the Breakpoint File, Local Assembly performed with `MindTheGap fill` and Post-Processing of the assembled sequences' file obtained.

    Args:
        - current_gap: str
            current gap/target identification
        - gfaFile: file
            GFA file containing the gaps' coordinates
        - chunkSize: int
            size of the chunk/flank region
        - extSize: int
            size of the gap/target extension on both sides (bp); determine start/end of the local assembly
        - maxLength: int
            maximum assembly length (bp)
        - minLength : int
            minimum assembly length (bp)
        - kmerSizeList: list
            list of k-mer sizes values
        - abundanceThresholdList: list
            list of minimal abundance thresholds values for solid k-mers
        - maxNodes: int
            maximum number of nodes in contig graph
        - nbCores: int
            number of cores 
        - maxMemory: int
            maximum memory for graph building (in MBytes)
        - verbosity: int
            verbosity level

    Return:
        - gapfillingFile: file
            file containing the obtained assembled sequence(s)
    """
    #----------------------------------------------------
    # Pre-Processing
    #----------------------------------------------------
    try:
        # Go in the 'outDir' directory.
        try:
            os.chdir(main.outDir)
        except OSError as err:
            print("File 'localAssemblyDBG.py': Something wrong with specified directory 'outDir'. \nOSError-{}".format(err))
            sys.exit(1)

        # Open the input GFA file.
        gfa = gfapy.Gfa.from_file(gfaFile)
        if not gfa:
            print("File 'localAssemblyDBG.py', function 'localAssemblyWithDBGAlgorithm()': Unable to open the input GFA file {}.".format(str(gfaFile)), file=sys.stderr)
            sys.exit(1)
        
        # Get the corresponding Gap line ('G' line).
        for _gap_ in gfa.gaps:
            if str(_gap_) == current_gap:
                current_gap = _gap_
                ## Create the object 'gap' from the class 'Gap'
                gap = Gap(current_gap)
                if not gap:
                    print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to create the object 'gap' from the class 'Gap'.", file=sys.stderr)
                    sys.exit(1)

        # Get some information on the current gap/target we are working on.
        gapLabel = gap.label()

        # Create two objects ('leftScaffold' and 'rightScaffold') from the class 'Scaffold'.
        leftScaffold = Scaffold(current_gap, gap.left, gfaFile)
        if not leftScaffold:
            print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to create the object 'leftScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)
        rightScaffold = Scaffold(current_gap, gap.right, gfaFile)
        if not rightScaffold:
            print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to create the object 'rightScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)

        # Get the gap/target flanking sequences (e.g. the flanking contigs sequences).
        leftFlankingSeq = str(leftScaffold.sequence())
        if not leftFlankingSeq:
            print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to get the left flanking sequence.", file=sys.stderr)
            sys.exit(1)
        rightFlankingSeq = str(rightScaffold.sequence())
        if not rightFlankingSeq:
            print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to get the right flanking sequence.", file=sys.stderr)
            sys.exit(1)

        # If chunk/flank size larger than length of scaffold(s), set the chunk/flank size to the minimal scaffold length.
        ## Left chunk/flank
        if chunkSize > leftScaffold.slen:
            chunk_L = leftScaffold.slen
        else:
            chunk_L = chunkSize
        ## Right chunk/flank
        if chunkSize > rightScaffold.slen:
            chunk_R = rightScaffold.slen
        else:
            chunk_R = chunkSize

        # Obtain the left and right regions for further getting the most represented flanking k-mers.
        ## Left region
        leftRegion = leftScaffold.chunk(chunk_L)
        if not leftRegion:
            print("File 'localAssemblyDBG.py': Unable to obtain the left region (left flank).", file=sys.stderr)
            sys.exit(1)
        ## Right region
        rightRegion = rightScaffold.chunk(chunk_R)
        if not rightRegion:
            print("File 'localAssemblyDBG.py': Unable to obtain the right region (right flank).", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print("File 'localAssemblyDBG.py': Something wrong with the 'Pre-Processing' step of the function 'localAssemblyWithDBGAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Breakpoint File (with offset of size k removed)
    #----------------------------------------------------
    # Go in the 'assemblyDir' directory.
    try:
        os.chdir(main.assemblyDir)
    except OSError as err:
        print("File 'localAssemblyDBG.py': Something wrong with specified directory 'assemblyDir'. \nOSError-{}".format(err))
        sys.exit(1)

    # Iterate over the k-mer sizes values, starting with the highest one.
    ## The list provided by the user must be a decreasing list.
    if len(kmerSizeList) == 0:
        print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Error with the input 'kmerSizeList' {}: empty list.".format(str(kmerSizeList)), file=sys.stderr)
        sys.exit(1)
    for k in kmerSizeList:

        try:
            # Left kmer.
            if leftScaffold.orient == "+":
                leftKmerRegion_start = int(str(leftRegion).split('-')[-1]) - extSize - k
                leftKmerRegion_end = int(str(leftRegion).split('-')[-1]) - extSize

            if leftScaffold.orient == "-":
                leftKmerRegion_start = int(str(leftRegion).split(':')[1].split('-')[0]) + extSize
                leftKmerRegion_end = int(str(leftRegion).split(':')[1].split('-')[0]) + extSize + k
                
            leftKmerRegion = str(leftRegion).split(':')[0] +":"+ str(leftKmerRegion_start) +"-"+ str(leftKmerRegion_end)
            leftKmer = getMostRepresentedKmer(main.bamFile, leftKmerRegion, leftScaffold.orient, k)
            if not leftKmer:
                print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to get the left kmer for {}.".format(str(leftRegion).split(':')[0]), file=sys.stderr)

            # Right kmer.
            if rightScaffold.orient == "+":
                rightKmerRegion_start = int(str(rightRegion).split(':')[1].split('-')[0]) + extSize
                rightKmerRegion_end = int(str(rightRegion).split(':')[1].split('-')[0]) + extSize + k

            if rightScaffold.orient == "-":
                rightKmerRegion_start = int(str(rightRegion).split('-')[-1]) - extSize - k
                rightKmerRegion_end = int(str(rightRegion).split('-')[-1]) - extSize
                
            rightKmerRegion = str(rightRegion).split(':')[0] +":"+ str(rightKmerRegion_start) +"-"+ str(rightKmerRegion_end)
            rightKmer = getMostRepresentedKmer(main.bamFile, rightKmerRegion, rightScaffold.orient, k)
            if not rightKmer:
                print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to get the right kmer for {}.".format(str(rightRegion).split(':')[0]), file=sys.stderr)

            # Reverse Left kmer.
            if rightScaffold.orient == "+":
                revLeftKmer = str(rc(rightFlankingSeq)[(len(rightFlankingSeq) - extSize - k):(len(rightFlankingSeq) - extSize)])
            if rightScaffold.orient == "-":
                revLeftKmer = str(rightFlankingSeq[(len(rightFlankingSeq) - extSize - k):(len(rightFlankingSeq) - extSize)])

            # Reverse Right kmer.
            if leftScaffold.orient == "+":
                revRightKmer = str(rc(leftFlankingSeq)[extSize:(extSize + k)])
            if leftScaffold.orient == "-":
                revRightKmer = str(leftFlankingSeq[extSize:(extSize + k)])

            # Get a breakpoint file containing the input sequences for the local assembly with `MindTheGap fill` (start and stop kmers).
            gfa_name = gfaFile.split('/')[-1]
            bkptFile = "{}.{}.g{}.flank{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, str(gapLabel), gap.length, chunkSize, k)
            try:
                with open(bkptFile, "w") as bkpt:
                    # Left kmer and Reverse Right kmer (dependent on the orientation of the left scaffold).
                    line1 = ">bkpt1_TargetID.{}_TargetLen.{} left_kmer.{}_len.{} offset_rm\n".format(str(gapLabel), gap.length, leftScaffold.name, k)
                    line2 = leftKmer
                    line7 = "\n>bkpt2_TargetID.{}_TargetLen.{} right_kmer.{}_len.{} offset_rm\n".format(str(gapLabel), gap.length, leftScaffold.name, k)
                    line8 = revRightKmer
                    
                    # Right kmer and Reverse Left kmer (dependent on the orientation of the right scaffold).
                    line3 = "\n>bkpt1_TargetID.{}_TargetLen.{} right_kmer.{}_len.{} offset_rm\n".format(str(gapLabel), gap.length, rightScaffold.name, k)
                    line4 = rightKmer
                    line5 = "\n>bkpt2_TargetID.{}_TargetLen.{} left_kmer.{}_len.{} offset_rm\n".format(str(gapLabel), gap.length, rightScaffold.name, k)
                    line6 = revLeftKmer
                    
                    bkpt.writelines([line1, line2, line3, line4, line5, line6, line7, line8])

            except IOError as err:
                print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Unable to open or write to the 'bkptFile' {}. \nIOError-{}".format(bkptFile, err))
                sys.exit(1)

        except Exception as e:
            print("File 'localAssemblyDBG.py': Something wrong with the 'Breakpoint File' creation step of the function 'localAssemblyWithDBGAlgorithm()'")
            print("Exception-{}".format(e))
            sys.exit(1)


        #----------------------------------------------------
        # Local Assembly performed with `MindTheGap fill`
        #----------------------------------------------------
        # Iterate over the abundance threshold values, starting with the highest one.
        ## The list provided by the user must be a decreasing list.
        if len(abundanceThresholdList) == 0:
            print("File 'localAssemblyDBG.py, function 'localAssemblyWithDBGAlgorithm()': Error with the input 'abundanceThresholdList' {}: empty list.".format(str(abundanceThresholdList)), file=sys.stderr)
            sys.exit(1)
        for a in abundanceThresholdList:

            try:
                print("LOCAL ASSEMBLY OF: {} for k={} and a={} (union)".format(str(gapLabel), k, a))
                
                # Input reads file containing the subsample of reads extracted during the 'Read Subsampling' step.
                try:
                    unionReadsFile = "{}.{}.g{}.flank{}.occ{}.bxu.fastq".format(gfa_name, str(gapLabel), gap.length, chunkSize, main.barcodesMinOcc)
                    subReadsFile = os.path.join(main.subsamplingDir, unionReadsFile)
                except FileNotFoundError as err:
                    print("File 'localAssemblyDBG.py': The input reads file {} doesn't exist. \nFileNotFoundError-{}".format(subReadsFile, err))
                    sys.exit(1)
                
                # Prefix of the output files containing the results obtained from `MindTheGap fill`. 
                outputPrefix = "{}.{}.g{}.flank{}.occ{}.k{}.a{}.bxu".format(gfa_name, str(gapLabel), gap.length, chunkSize, main.barcodesMinOcc, k, a)

                # Determine the maximum assembly length (bp) if the gap/target length is known.
                ## NB: if the gap/target length is unknown, it is set to 0.
                ## Add twice the extension size to the gap/target length (extension on both sides of the gap/target) and twice the length of reads (e.g. 2x 150 bp) to be large
                if gap.length >= maxLength:
                    maxLength = gap.length + 2*extSize + 2*150

                # Perform the local assembly step with `MindTheGap fill`.
                res, success = fillGapWithMTGFill(gapLabel, subReadsFile, bkptFile, k, a, maxLength, maxNodes, nbCores, maxMemory, verbosity, outputPrefix)

                # For each assembled sequence, check that the assembly length is larger than the minimum assembly length required, and if so, add the assembly to the 'assembliesFile'.
                if success:
                    success = False

                    # Save the sequences having an assembly length larger than the minimum assembly length required in the file 'assembliesFile'.
                    try:
                        try:
                            assembliesFile = os.path.abspath(main.assemblyDir +"/"+ outputPrefix + ".insertions_filtered.fasta")
                        except FileNotFoundError as err:
                            print("File 'localAssemblyDBG.py': The output 'assembliesFile' {} doesn't exist. \nFileNotFoundError-{}".format(str(assembliesFile), err))
                            sys.exit(1)
                        with open(res, "r") as mtgOutput, open(assembliesFile, "w") as assemblies:
                            records = SeqIO.parse(mtgOutput, "fasta")
                            for record in records:
                                if len(record.seq) >= minLength:
                                    SeqIO.write(record, assemblies, "fasta")
                                    success = True
                    except IOError as err:
                        print("File 'localAssemblyDBG.py', function 'localAssemblyWithDBGAlgorithm()': Unable to open the file {} or to open/write to the file {}. \nIOError-{}".format(str(res), str(assembliesFile), err))
                        sys.exit(1)
                    if not success:
                        res = "Local Assembly not completed for k{} and lower...".format(str(max(kmerSizeList)))

            except Exception as e:
                print("File 'localAssemblyDBG.py': Something wrong with the 'Local Assembly performed with `MindTheGap fill`' step of the function 'localAssemblyWithDBGAlgorithm()'")
                print("Exception-{}".format(e))
                sys.exit(1)
            
            # Case of successful local assembly, break the loop on 'a'.
            if success:
                break

        # Case of successful local assembly, break the loop on 'k'.
        if success:
            break


    #----------------------------------------------------
    # Post-Processing
    #----------------------------------------------------
    try:
        # Output file containing the assembled sequence(s).
        insertionsFile = "{}.{}.g{}.flank{}.occ{}.k{}.a{}.bxu.insertions_filtered.fasta".format(gfa_name, str(gapLabel), gap.length, chunkSize, main.barcodesMinOcc, k, a)

        # Case of unsuccessful local assembly.
        if not success:
            print("\n{}: {}\n".format(gapLabel, res))
            try:
                gapfillingFile = os.path.abspath(insertionsFile)
                open(gapfillingFile, 'w').close()
            except FileNotFoundError as err:
                print("File 'localAssemblyDBG.py': The output 'gapfillingFile' {} doesn't exist. \nFileNotFoundError-{}".format(str(gapfillingFile), err))
                sys.exit(1)

        # Case of successful local assembly.
        if success:
            print("\n{}: Successful Local Assembly for k{} !\n". format(gapLabel, str(k)))

            # Save and pre-process the file containing the assembled sequence(s) for further qualitative evaluation.
            ## Modify the 'insertionFile' and save it to a new file ('gapfillingFile') so that the 'solution x/y' part appears in record.id (and not just in record.description)
            try:
                gapfillingFile = os.path.abspath(main.assemblyDir +"/"+ outputPrefix + "..insertions_filtered.fasta")
            except FileNotFoundError as err:
                print("File 'localAssemblyDBG.py': The output 'gapfillingFile' {} doesn't exist. \nFileNotFoundError-{}".format(str(gapfillingFile), err))
                sys.exit(1)
            try:
                with open(insertionsFile, "r") as original, open(gapfillingFile, "w") as corrected:
                    records = SeqIO.parse(original, "fasta")
                    for record in records:
                        if "solution" in record.description:
                            record.id = record.id + "_sol_" + record.description.split(" ")[-1]
                        else:
                            record.id = record.id + "_sol_1/1"
                        SeqIO.write(record, corrected, "fasta")
            except IOError as err:
                print("File 'localAssemblyDBG.py', function 'localAssemblyWithDBGAlgorithm()': Unable to open the file {} or to open/write to the file {}. \nIOError-{}".format(str(insertionsFile), str(gapfillingFile), err))
                sys.exit(1)
    
    except Exception as e:
        print("File 'localAssemblyDBG.py': Something wrong with the 'Post-Processing' step of the function 'localAssemblyWithDBGAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    return gapfillingFile

