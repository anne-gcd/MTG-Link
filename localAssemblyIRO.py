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

"""Module 'localAssemblyIRO.py': Local Assembly With the IRO (Iterative Read Overlap) algorithm

The module 'localAssemblyIRO.py' enables to perform the local assembly step of the MTG-Link gap-filling pipeline, using an Iterative Read Overlap (IRO) algorithm. 
The IRO algorithm is based on-the-fly computations of read overlaps and iterative extensions of the current assembly sequence, using the subsample of linked reads obtained during the 'Read Subsampling' step of the MTG-Link pipeline.
"""

from __future__ import print_function
import collections
import itertools
import os
import regex
import sys
import gfapy
import time
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter, pos
import main
from helpers import Gap, Scaffold
from ProgDynOptim import DynamicMatrixOptim


# Increase the maximum recursion depth in Python.
sys.setrecursionlimit(50000)


#----------------------------------------------------
# indexReadBySeed function
#----------------------------------------------------
def indexReadBySeed(read, pos, rcRead, seedSize, seedDict):
    """
    To index a read by its seed.

    Args:
        - read: str
            sequence of the current read to index
        - pos: str
            position of the current read in readsList (list containing all reads' sequences)
        - rcRead: str
            sequence of the reverse complement of the current read
        - seedSize: int
            size of the seed used for indexing the reads
        - seedDict: dict
            this function will output a dictionary 'seedDict' containing the reads indexed by their seed
            key = seed's sequence ; value = list of positions of reads having this seed in readList

    Output:
        - seedDict: dict
            dictionary 'seedDict' containing the reads indexed by their seed
            key = seed's sequence ; value = list of positions of reads having this seed in readList
    """
    try:
        # Index read by its seed.
        seed = read[:seedSize]
        if seed in seedDict:
            seedDict[seed].append(str(pos))
        else:
            seedDict[seed] = [str(pos)]

        # Index reverse complement of read by its seed as well.
        seed = rcRead[:seedSize]
        if seed in seedDict:
            seedDict[seed].append("-"+str(pos))
        else:
            seedDict[seed] = ["-"+str(pos)]
        
        if len(seedDict) == 0:
            print("File 'localAssemblyIRO.py, function 'indexReadBySeed()': Error with the output 'seedDict' {}: empty dict.".format(str(seedDict)), file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the function 'indexReadBySeed()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# findOverlappingReads function
#----------------------------------------------------
def findOverlappingReads(assembly, lenRead, readsList, seedSize, minOverlapSize, dmax, seedDict):
    """
    To find the reads overlapping with the current assembly's sequence S.

    Args:
        - assembly: str
            current assembly's sequence
        - lenRead: int
            length of the read from which we want to extend
        - readsList: list
            reads list used for the local assembly step with the IRO algorithm
        - seedSize: int
            size of the seed used for indexing the reads
        - minOverlapSize: int
            minimum overlapping size
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads
        - seedDict: dict
            dictionary 'seedDict' containing reads indexed by their seed
            key = seed's sequence ; value = list of positions of reads having this seed in readList

    Returns:
        - overlappingReads: list
            list 'overlappingReads' containing all the overlapping reads' sequences, along with the index of the beginning of the overlap on S, and the index of the beginning of extension on R,
            referenced as [read's sequence, index of beginning of overlap, index of beginning of extension]
            this list is sorted automatically by smallest 'posOverlap', e.g. by largest overlap
    """
    try:
        overlappingReads = []

        # Get the putative reads (e.g. reads having a seed onto the current assembly's sequence).
        for posOverlap in range(len(assembly) - lenRead + 1, len(assembly) - minOverlapSize - seedSize):
            seed = assembly[posOverlap : posOverlap+seedSize]
            if not seed:
                print("File 'localAssemblyIRO.py, function 'findOverlappingReads()': Error with the current seed determination.", file=sys.stderr)
                sys.exit(1)

            if seed in seedDict:
                putativeReads = seedDict[seed]

                # For each putative read, search for an overlap between the current assembly's sequence and the putative read, using the dynamic programmation.
                for putRead in putativeReads:

                    # Get the sequence of the read.
                    if len(readsList) == 0:
                        print("File 'localAssemblyIRO.py, function 'findOverlappingReads()': Error with the input 'readsList' {}: empty list.".format(str(readsList)), file=sys.stderr)
                        sys.exit(1)
                    if '-' in str(putRead):
                        read = str(Seq(readsList[int(putRead.split('-')[1])]).reverse_complement())
                    else:
                        read = readsList[int(putRead)]
                    
                    # Perform the alignment with the optimized dynamic programmation.
                    dm = DynamicMatrixOptim(assembly[(posOverlap+seedSize):], read[seedSize:], dmax)
                    dist, posR = dm.getEditDistanceAndGenomePosition()
                    
                    # Overlap found.
                    if dist != None:
                        posExt = posR + seedSize
                        overlappingReads.append([read, posOverlap, posExt])

        return overlappingReads

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the function 'findOverlappingReads()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# extendReadWithOverlappingReads function
#----------------------------------------------------
def extendReadWithOverlappingReads(assembly, lenRead, inputSeqName, STOP, seedDict, assemblyHash, readsList, seedSize, minOverlapSize, abundanceMinList, dmax, maxLength, bpTotal, startTime, iroLog):
    """
    To extend a read's sequence with overlapping reads.
    This is an iterative function.
    The Boolean value it returns represents the success of the gap-filling
    NB: extGroup is a dictionary containing the extension's sequence as key, and the reads sharing this extension as value
        (value format: [read's sequence, index of beginning of overlap])

    Args:
        - assembly: str
            current assembly's sequence
        - lenRead: int
            length of the read from which we want to extend
        - inputSeqName: str
            name of the input k-mers sequences
        - STOP: str
            k-mer STOP (target)
        - seedDict: dict
            dictionary 'seedDict' containing reads indexed by their seed
            key = seed's sequence ; value = list of positions of reads having this seed in readList
        - assemblyHash = hashtable/dict
            hashtable/dictionary indicating if the search for overlapping reads has already been performed on the corresponding sequence (key):
            key = the last 70 bp of the current assembly's sequence ; value = Boolean value (0: overlapping reads search not performed / 1: overlapping reads search performed)
        - readsList: list
            reads list used for the local assembly step with the IRO algorithm
        - seedSize: int
            size of the seed used for indexing the reads
        - minOverlapSize: int
            minimum overlapping size
        - abundanceMinList: list
            list of minimal abundance(s) of reads used for gapfilling
            extension's groups having less than this number of reads are discarded from the graph
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads
        - maxLength: int
            maximum assembly length (bp)
        - bpTotal: int
            number of total bp added to all possible assembled sequences
        - startTime: time
            starting time of the exploration for extending the assembly sequence
        - iroLog: file
            the temporary (e.g. incomplete) gap-filled sequences will be saved in this log file

    Return:
        - str, Boolean
            - the gap-filled sequence (assembly) and a Boolean variable equal to True if a solution is found 
              (e.g. we arrived to STOP kmer (target), with at most 2 substitutions in the target sequence)
            OR
            - the reason why the gap-filling failed and a Boolean variable equal to False if no solution is found
        - bpTotal: the number of total bp added to all possible assembled sequences
    """
    try:
        # Base cases.

        ## Target reached (with at most 2 substitutions): Successful gap-filling.
        kmer_STOP = "({})".format(STOP)
        match = regex.findall(str(kmer_STOP)+"{s<=2}", assembly, overlapped=True)
        if match != []:
            return assembly, True, bpTotal
        
        ## Assembly length superior to 'maxLength' specified by the user: Gap-filling aborted.
        if len(assembly) > maxLength:
            return "|S| > maxLength", False, bpTotal

        ## Number of total bp added to all possible assembled sequences higher than 100*maxLength: Gap-filling aborted.
        if bpTotal > 10*maxLength:
            return "Too many explorations: No solution", False, bpTotal

        ## Exploration takes too much time (> 25% maxLength e.g. we allow 0.25 s per bp): Gap-filling aborted.
        if (time.time() - startTime) > 0.25*maxLength:
            return "Exploration takes too much time: No solution", False, bpTotal

        ## Path already explored.
        if len(assembly) >= 70:
            # Check that we didn't already search for overlapping reads on this region (e.g. on the last 70 bp of the current assembly's sequence).
            if assemblyHash[assembly[-70:]] == 1:
                return "Path already explored: No solution", False, bpTotal
                
        # Search for reads overlapping with the current assembly's sequence.
        if len(readsList) == 0:
            print("File 'localAssemblyIRO.py, function 'extendReadWithOverlappingReads()': Error with the input 'readsList' {}: empty list.".format(str(readsList)), file=sys.stderr)
            sys.exit(1)
        overlappingReads = findOverlappingReads(assembly, lenRead, readsList, seedSize, minOverlapSize, dmax, seedDict)
        
        if not overlappingReads:
            try:
                with open(iroLog, "a") as log:
                    log.write("\n>" + inputSeqName + " _ No.overlapping.reads")
                    log.write("\n"+str(assembly)+"\n")
            except IOError as err:
                print("File 'localAssemblyIRO.py', function 'extendReadWithOverlappingReads()': Unable to open or write to the file {}. \nIOError-{}".format(str(iroLog), err))
                sys.exit(1)
            return "No overlapping reads", False, bpTotal

        # Group the overlapping reads by their extension.
        extGroup = {}

        # Populate extGroup.
        '''NB: 'overlappingReads' list sorted automatically by smallest 'posOverlap', e.g. by largest overlap'''
        for (readSeq, posOverlap, posExt) in overlappingReads:

            # If no extension, don't add it to extGroup.
            if readSeq[posExt:] == "":
                continue

            # Add first extension to extGroup.
            if len(extGroup) == 0:
                extGroup[readSeq[posExt:]] = [[readSeq, posOverlap]]

            # Add all extensions to extGroup.
            elif len(extGroup) > 0:
                for extension in extGroup:

                    # Check that current extension is smaller than the one(s) in extGroup.
                    if len(readSeq[posExt:]) < len(extension):
                        # Current extension already in extGroup.
                        if readSeq[posExt:] == extension[:len(readSeq[posExt:])]:
                            new_extension = readSeq[posExt:]
                            extGroup[new_extension] = extGroup[extension]
                            extGroup[new_extension].append([readSeq, posOverlap])
                            del extGroup[extension]
                            added_to_extGroup = True
                            break
                        # Current extension is partially in extGroup.
                        elif readSeq[posExt] == extension[0]:
                            i = 1
                            while i < len(readSeq[posExt:]):
                                if readSeq[posExt+i] == extension[i]:
                                    i += 1
                                else:
                                    break
                            new_extension = extension[:i]
                            extGroup[new_extension] = extGroup[extension]
                            extGroup[new_extension].append([readSeq, posOverlap])
                            del extGroup[extension]
                            added_to_extGroup = True
                            break
                        # Current extension not already in extGroup.
                        else:
                            added_to_extGroup = False

                    # Current extension is not smaller than the one(s) in extGroup.
                    else:
                        # Current extension already in extGroup.
                        if readSeq[posExt:posExt+len(extension)] == extension:
                            extGroup[extension].append([readSeq, posOverlap])
                            added_to_extGroup = True
                            break
                        # Current extension is partially in extGroup.
                        elif readSeq[posExt] == extension[0]:
                            i = 1
                            while i < len(extension):
                                if readSeq[posExt+i] == extension[i]:
                                    i += 1
                                else:
                                    break
                            new_extension = extension[:i]
                            extGroup[new_extension] = extGroup[extension]
                            extGroup[new_extension].append([readSeq, posOverlap])
                            del extGroup[extension]
                            added_to_extGroup = True
                            break
                        # Current extension not already in extGroup.
                        else:
                            added_to_extGroup = False
                            
                # Current extension not already in extGroup.
                if not added_to_extGroup:
                    extGroup[readSeq[posExt:]] = [[readSeq, posOverlap]]

            # Sort extGroup by the smallest extension.
            extGroup = collections.OrderedDict(sorted(extGroup.items(), key=lambda t: len(t[0])))

        # Update 'assemblyHash' to indicate that we performed the search for overlapping reads on this region (e.g. on the last 70 bp of the current assembly's sequence).
        assemblyHash[assembly[-70:]] = 1

        # Filter extGroup by the number of reads sharing an extension (argument 'abundanceMin').
        if len(abundanceMinList) == 0:
            print("File 'localAssemblyIRO.py, function 'extendReadWithOverlappingReads()': Error with the input 'abundanceMinList' {}: empty list.".format(str(abundanceMinList)), file=sys.stderr)
            sys.exit(1)
        for abundanceMin in abundanceMinList:
            extGroupFiltered = extGroup.copy()
            if not extGroupFiltered:
                print("File 'localAssemblyIRO.py, function 'extendReadWithOverlappingReads()': Error with copying 'extGroup' to 'extGroupFiltered'.", file=sys.stderr)
                sys.exit(1)

            for extension in list(extGroupFiltered.keys()):
                if len(extGroupFiltered[extension]) < abundanceMin:
                    del extGroupFiltered[extension]
            
            # Iterate over the values of abundanceMin only if the number of reads sharing an extension < 'abundanceMin'.
            if not extGroupFiltered:
                continue
            else:
                break

        # If number of reads sharing an extension < minimal 'abundanceMin' provided, stop the extension.
        if not extGroupFiltered:
            try:
                with open(iroLog, "a") as log:
                    log.write("\n>" + inputSeqName + " _ No.extension")
                    log.write("\n"+str(assembly)+"\n")
            except IOError as err:
                print("File 'localAssemblyIRO.py', function 'extendReadWithOverlappingReads()': Unable to open or write to the file {}. \nIOError-{}".format(str(iroLog), err))
                sys.exit(1)
            return "No extension", False, bpTotal

        # Sort extGroup by the extension whose read has the largest overlap with the current assembly's sequence (smallest 'posOverlap'). 
        '''NB: values of extGroup sorted by reads having the larger overlap'''
        extGroupFiltered = collections.OrderedDict(sorted(extGroupFiltered.items(), key=lambda  t: t[1][0][1]))

        # If there are many possible extensions, iterate over the 5 first extensions.
        if len(extGroupFiltered) > 5:
            finalExtGroup = dict(itertools.islice(extGroupFiltered.items(), 5))
        else:
            finalExtGroup = extGroupFiltered.copy()
            if not finalExtGroup:
                print("File 'localAssemblyIRO.py, function 'extendReadWithOverlappingReads()': Error with copying 'extGroupFiltered' to 'finalExtGroup'.", file=sys.stderr)
                sys.exit(1)

        # Iterative extension of the assembly's sequence S.
        for extension in finalExtGroup:

            # Update 'assemblyHash' with the new region for which we will search for overlapping reads (with value '0' if search not already performed, or with value '1' if search already performed).
            if (assembly+extension)[-70:] in assemblyHash.keys():
                assemblyHash[(assembly+extension)[-70:]] = 1
            else:
                assemblyHash[(assembly+extension)[-70:]] = 0
            
            # Extend the updated assembly sequence (e.g. assembly+extension sequence).
            bpTotal = bpTotal + len(extension)
            res, success, bpTotal = extendReadWithOverlappingReads(assembly+extension, len(extGroupFiltered[extension][0][0]), inputSeqName, STOP, seedDict, assemblyHash, readsList, seedSize, minOverlapSize, abundanceMinList, dmax, maxLength, bpTotal, startTime, iroLog)

            # If we find a complete gap-filled sequence (e.g. we reach the kmer STOP), return the assembly sequence along with True.
            if success:
                return res, True, bpTotal

        # If we don't find a complete gap-filled sequence, return the reason why the gap-filling was not successful along with False.
        return res, False, bpTotal

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the function 'extendReadWithOverlappingReads()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# fillGapWithIROAlgo function
#----------------------------------------------------
def fillGapWithIROAlgo(gapLabel, readsList, bkptFile, seedSize, minOverlapSize, abundanceMinList, dmax, maxLength):
    """
    To execute the IRO algorithm, that is based on on-the-fly computations of read overlaps and iterative extensions of the current assembly sequence.
    The local assembly step is performed between the sequences surrounding the extended gap (e.g. the kmers of the breakpoint file), in the forward orientation.

    Args:
        - gapLabel: str
            label of the gap
        - readsList: list
            reads list used for the local assembly step with the IRO algorithm
        - bkptFile: file
            breakpoint file containing the input sequences for the local assembly with the IRO algorithm (START and STOP kmers)
            (breakpoint sequences with offset of size ext removed)
        - seedSize: int
            size of the seed used for indexing the reads
        - minOverlapSize: int
            minimum overlapping size
        - abundanceMinList: list
            list of minimal abundance(s) of reads used for gapfilling
            extension's groups having less than this number of reads are discarded from the graph
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads
        - maxLength: int
            maximum assembly length (bp)

    Return:
        - res: str
            - the assembled sequence from the read containing the k-mer START to the read containing the k-mer STOP (assembly) (and so the Boolean value equal to True)
              if a solution is found (e.g. we arrived to STOP kmer (target), with at most 2 substitutions in the target sequence)
            OR
            - the reason why the gap-filling failed (and so the Boolean value equal to False)
              if no solution is found
        - success: boolean
            True/False
    """
    try:
        # Initiate the four main variables.
        seedDict = {}
        readsWithStartList = []
        positionReadInReadsList = 0
        assemblyHash = {}

        # Initiate the timeout.
        bpTotal = 0
        startTime = time.time()

        # Initiate the log file.
        iroLog = str(gapLabel) + "_IROAlgo.log"

        # Get the k-mers gap flanking sequences (kmers START and STOP) (e.g. source and target).
        try:
            with open(bkptFile, "r") as bkpt:
                for record in SeqIO.parse(bkpt, "fasta"):
                    if (record.id == "start") or (record.id == "START") or ("left" in record.description):
                        START = str(record.seq)
                        inputSeqName = record.id
                    if (record.id == "stop") or (record.id == "STOP") or ("right" in record.description):
                        STOP = str(record.seq)
                        if record.id != inputSeqName:
                            inputSeqName += "-" + record.id
        except IOError as err:
            print("File 'localAssemblyIRO.py', function 'fillGapWithIROAlgo()': Unable to open the file {}. \nIOError-{}".format(str(bkptFile), err))
            sys.exit(1)

        if not START:
            print("File 'localAssemblyIRO.py, function 'fillGapWithIROAlgo()': Error with the creation of the START k-mer.", file=sys.stderr)
            sys.exit(1)
        if not STOP:
            print("File 'localAssemblyIRO.py, function 'fillGapWithIROAlgo()': Error with the creation of the STOP k-mer.", file=sys.stderr)
            sys.exit(1)

        # Iterate over the reads of 'readsList' to obtain the 'seedDict' dictionary and the 'readsWithStartList' list.
        if len(readsList) == 0:
            print("File 'localAssemblyIRO.py, function 'fillGapWithIROAlgo()': Error with the input 'readsList' {}: empty list.".format(str(readsList)), file=sys.stderr)
            sys.exit(1)
        for read in readsList:

            # Get the reverse complement of the current read.
            rcRead = str(Seq(read).reverse_complement())
            # Seed the read and update the 'seedDict' dictionary.
            indexReadBySeed(read, positionReadInReadsList, rcRead, seedSize, seedDict)

            # Search if the read contains the whole kmer START's sequence (e.g. source k-mer) and update the 'readsWithStartList' list.
            if START in read:
                readsWithStartList.append([str(positionReadInReadsList), read.index(START)])
            elif START in rcRead:
                readsWithStartList.append(["-"+str(positionReadInReadsList), rcRead.index(START)])

            # Increment the position of the current read in 'readList'
            positionReadInReadsList += 1

        # Sort the 'readsWithStartList' list by the minimum extension size (e.g. by the maximum index).
        readsWithStartList = sorted(readsWithStartList, key=itemgetter(1), reverse=True)
        
        # If there is no read containing the kmer START, raise an exception.
        if not readsWithStartList:
            return "No read in the dataset provided contains the kmer START", False

        # Extend the reads containing the whole kmer START's sequence.
        for (posReadInReadsList, index) in readsWithStartList:

            # Get the sequence of the read.
            if '-' in str(posReadInReadsList):
                read = str(Seq(readsList[int(posReadInReadsList.split('-')[1])]).reverse_complement())
            else:
                read = readsList[int(posReadInReadsList)]

            # Extend the assembly sequence (e.g. the current read containing the whole kmer START's sequence).
            assemblyHash[read[-70:]] = 0
            bpTotal = bpTotal + len(read)
            if len(abundanceMinList) == 0:
                print("File 'localAssemblyIRO.py, function 'fillGapWithIROAlgo()': Error with the input 'abundanceMinList' {}: empty list.".format(str(abundanceMinList)), file=sys.stderr)
                sys.exit(1)
            res, success, bpTotal = extendReadWithOverlappingReads(read, len(read), inputSeqName, STOP, seedDict, assemblyHash, readsList, seedSize, minOverlapSize, abundanceMinList, dmax, maxLength, bpTotal, startTime, iroLog)

            # Case of successful gap-filling.
            break

        return res, success

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the function 'fillGapWithIROAlgo()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# localAssemblyWithIROAlgorithm function
#----------------------------------------------------
def localAssemblyWithIROAlgorithm(current_gap, gfaFile, chunkSize, extSize, maxLength, seedSize, minOverlapSize, abundanceMinList, dmax):
    """
    To perform the Local Assembly step using a IRO (Iterative Read Overlap) algorithm.  
    The IRO algorithm is based on on-the-fly computations of read overlaps and iterative extensions of the current assembly sequence. This module is executed on the subsample of reads retrieved during the 'Read Subsampling' step.
    This consists of five main steps: Pre-processing of the current gap, getting the Breakpoint File, getting the 'readList' list, Local Assembly performed with the IRO algorithm and Post-Processing of the gap-filled sequences obtained.

    Args:
        - current_gap: str
            current gap identification
        - gfaFile: file
            GFA file containing the gaps' coordinates
        - chunkSize: int
            size of the chunk region
        - extSize: int
            size of the gap extension on both sides (bp); determine start/end of the local assembly
        - maxLength: int
            maximum assembly length (bp)
        - seedSize: int
            size of the seed used for indexing the reads
        - minOverlapSize: int
            minimum overlapping size for reads overlaps
        - abundanceMinList: list
            list of minimal abundance(s) of reads used for gapfilling
            extension's groups having less than this number of reads are discarded from the graph
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads

    Return:
        - gapfillingFile: file
            file containing the obtained gap-filled sequence
    """
    #----------------------------------------------------
    # Pre-Processing
    #----------------------------------------------------
    try:
        # Go in the 'outDir' directory.
        try:
            os.chdir(main.outDir)
        except OSError as err:
            print("File 'localAssemblyIRO.py': Something wrong with specified directory 'outDir'. \nOSError-{}".format(err))
            sys.exit(1)

        # Open the input GFA file.
        gfa = gfapy.Gfa.from_file(gfaFile)
        if not gfa:
            print("File 'localAssemblyIRO.py', function 'localAssemblyWithIROAlgorithm()': Unable to open the input GFA file {}.".format(str(gfaFile)), file=sys.stderr)
            sys.exit(1)
        
        # Get the corresponding Gap line ('G' line).
        for _gap_ in gfa.gaps:
            if str(_gap_) == current_gap:
                current_gap = _gap_
                ## Create the object 'gap' from the class 'Gap'
                gap = Gap(current_gap)
                if not gap:
                    print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Unable to create the object 'gap' from the class 'Gap'.", file=sys.stderr)
                    sys.exit(1)

        # Get some information on the current gap we are working on.
        gapLabel = gap.label()

        # Create two objects ('leftScaffold' and 'rightScaffold') from the class 'Scaffold'.
        leftScaffold = Scaffold(current_gap, gap.left, gfaFile)
        if not leftScaffold:
            print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Unable to create the object 'leftScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)
        rightScaffold = Scaffold(current_gap, gap.right, gfaFile)
        if not rightScaffold:
            print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Unable to create the object 'rightScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)

        # Get the gap flanking sequences (e.g. the flanking contigs sequences).
        leftFlankingSeq = str(leftScaffold.sequence())
        if not leftFlankingSeq:
            print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Unable to get the left flanking sequence.", file=sys.stderr)
            sys.exit(1)
        rightFlankingSeq = str(rightScaffold.sequence())
        if not rightFlankingSeq:
            print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Unable to get the right flanking sequence.", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the 'Pre-Processing' step of the function 'localAssemblyWithIROAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Breakpoint File (with offset of size ext removed)
    #----------------------------------------------------
    # Go in the 'assemblyDir' directory.
    try:
        os.chdir(main.assemblyDir)
    except OSError as err:
        print("File 'localAssemblyIRO.py': Something wrong with specified directory 'assemblyDir'. \nOSError-{}".format(err))
        sys.exit(1)

    try:
        # Get a breakpoint file containing the input sequences for the local assembly with the IRO algorithm (start and stop kmers).
        gfa_name = gfaFile.split('/')[-1]
        bkptFile = "{}.{}.g{}.c{}.kmersStartStop.fasta".format(gfa_name, str(gapLabel), gap.length, chunkSize)
        try:
            with open(bkptFile, "w") as bkpt:

                # Start sequence (left kmer).
                line1 = ">ctg{}_START left_kmer_len.31\n".format(leftScaffold.scaffold)
                line2 = leftFlankingSeq[(leftScaffold.slen - extSize - 31):(leftScaffold.slen - extSize)]
                
                # Stop sequence (right kmer).
                line3 = "\n>ctg{}_STOP right_kmer_len.31\n".format(rightScaffold.scaffold)
                line4 = rightFlankingSeq[extSize:(extSize + 31)]

                bkpt.writelines([line1, line2, line3, line4])

        except IOError as err:
                print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Unable to open or write to the 'bkptFile' {}. \nIOError-{}".format(bkptFile, err))
                sys.exit(1)

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the 'Breakpoint File' creation step of the function 'localAssemblyWithIROAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # 'readsList' list
    #----------------------------------------------------
    try:
        # Input reads file containing the subsample of reads extracted during the 'Read Subsampling' step.
        try:
            unionReadsFile = "{}.{}.g{}.c{}.f{}.bxu.fastq".format(gfa_name, str(gapLabel), gap.length, chunkSize, main.barcodesMinFreq)
            subReadsFile = os.path.join(main.subsamplingDir, unionReadsFile)
        except FileNotFoundError as err:
            print("File 'localAssemblyIRO.py': The input reads file {} doesn't exist. \nFileNotFoundError-{}".format(str(subReadsFile), err))
            sys.exit(1)

        # Create the list 'readsList' containing all reads' sequences obtained from the 'Read Subsampling' step.
        readsList = []
        try:
            with open(subReadsFile, "r") as readsFile:
                readsList = [str(read.seq) for read in SeqIO.parse(subReadsFile, "fastq")]
        except IOError as err:
            print("File 'localAssemblyIRO.py', function 'localAssemblyWithIROAlgorithm()': Unable to open the input 'subReadsFile' used for the local assembly step {}. \nIOError-{}".format(str(subReadsFile), err))
            sys.exit(1)

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the 'readsList list' creation step of the function 'localAssemblyWithIROAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Local Assembly performed with the IRO algorithm
    #----------------------------------------------------
    try:
        # Convert the list of abundance min to a string.
        if len(abundanceMinList) == 0:
            print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Error with the input 'abundanceMinList' {}: empty list.".format(str(abundanceMinList)), file=sys.stderr)
            sys.exit(1)
        abundanceMinString = '-'.join(map(str, abundanceMinList))
        
        print("GAP-FILLING OF: {} for s={}, o={}, a={} and dmax={} (union)".format(str(gapLabel), seedSize, minOverlapSize, abundanceMinString, dmax))

        # Determine the maximum assembly length (bp) if the gap length is known.
        ## NB: if the gap length is unknown, it is set to 0.
        ## Add twice the extension size to the gap length (extension on both sides of the gap) and twice the length of reads (e.g. 2x 150 bp) to be large
        if gap.length >= maxLength:
            maxLength = gap.length + 2*extSize + 2*150

        # Perform the local assembly step with the IRO algorithm.
        res, success = fillGapWithIROAlgo(gapLabel, readsList, bkptFile, seedSize, minOverlapSize, abundanceMinList, dmax, maxLength)

    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the 'Local Assembly performed with the IRO algorithm' step of the function 'localAssemblyWithIROAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Post-Processing
    #----------------------------------------------------
    try:
        # Output file containing the gap-filled sequence(s).
        insertionsFile = "{}.{}.g{}.c{}.f{}.s{}.o{}.a{}.dmax{}.bxu.insertions.fasta".format(gfa_name, str(gapLabel), gap.length, chunkSize, main.barcodesMinFreq, seedSize, minOverlapSize, abundanceMinString, dmax)

        # Case of unsuccessful gap-filling.
        if not success:
            print("\n{}: {}\n".format(gapLabel, res))
            try:
                gapfillingFile = os.path.abspath(insertionsFile)
            except FileNotFoundError as err:
                print("File 'localAssemblyIRO.py': The output 'gapfillingFile' {} doesn't exist. \nFileNotFoundError-{}".format(str(gapfillingFile), err))
                sys.exit(1)
            # Initiate the gapfillingFile.
            ## Write nothing to it as there is no gap-filled sequence found.
            try:
                with open(gapfillingFile, "w") as f:
                    pass
            except IOError as err:
                print("File 'localAssemblyIRO.py', function 'localAssemblyWithIROAlgorithm()': Unable to open or write to the input 'gapfillingFile' {}. \nIOError-{}".format(str(gapfillingFile), err))
                sys.exit(1)
            
            # Save the reason why there is no complete gap-filling in a log file.
            iroLog = str(gapLabel) + "_IROresults.log"
            try:
                with open(iroLog, "a") as log:
                    log.write("\n" + res)
            except IOError as err:
                print("File 'localAssemblyIRO.py', function 'localAssemblyWithIROAlgorithm()': Unable to open or write to the file {}. \nIOError-{}".format(str(iroLog), err))
                sys.exit(1)

        # Case of successful gap-filling.
        if success:
            print("\n{}: Successful Gap-filling !\n". format(gapLabel))

            # Get the k-mers gap flanking sequences (kmers START and STOP).
            START = leftFlankingSeq[(leftScaffold.slen - extSize - 31):(leftScaffold.slen - extSize)]
            STOP = rightFlankingSeq[extSize:(extSize + 31)]
            inputSeqName = "ctg{}_START-ctg{}_STOP".format(gap.left, gap.right)

            # Get the target sequence found in the assembly, with at most 2 substitutions.
            kmer_STOP = "({})".format(STOP)
            match = regex.findall(str(kmer_STOP)+"{s<=2}", res, overlapped=True)
            STOP = match[0]
            
            if not START:
                print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Error with the creation of the final START k-mer.", file=sys.stderr)
                sys.exit(1)
            if not STOP:
                print("File 'localAssemblyIRO.py, function 'localAssemblyWithIROAlgorithm()': Error with the creation of the final STOP k-mer.", file=sys.stderr)
                sys.exit(1)
            
            # Save and pre-process the gap-filled sequence obtained for further qualitative evaluation.
            try:
                gapfillingFile = os.path.abspath(insertionsFile)
            except FileNotFoundError as err:
                print("File 'localAssemblyIRO.py': The output 'gapfillingFile' {} doesn't exist. \nFileNotFoundError-{}".format(str(gapfillingFile), err))
                sys.exit(1)
            try:
                with open(gapfillingFile, "a") as finalFile:
                    assemblyStartBegin = res.index(START)
                    assemblyStopBegin = res.index(STOP)
                    seq = res[assemblyStartBegin+len(START) : assemblyStopBegin]
                    seqName = "assembly." + inputSeqName + "_len_" + str(len(seq))
                    finalFile.write(">" + seqName)
                    finalFile.write("\n" + seq + "\n")
            except IOError as err:
                print("File 'localAssemblyIRO.py', function 'localAssemblyWithIROAlgorithm()': Unable to open or write to the file {}. \nIOError-{}".format(str(gapfillingFile), err))
                sys.exit(1)
    
    except Exception as e:
        print("File 'localAssemblyIRO.py': Something wrong with the 'Post-Processing' step of the function 'localAssemblyWithIROAlgorithm()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    return gapfillingFile
    
    