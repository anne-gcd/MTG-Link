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

"""Module 'IRO.py': Iterative Read Overlap (IRO) algorithm

The module 'IRO.py' enables to perform the local assembly step of the MTG-Link gap-filling pipeline, using an Iterative Read Overlap (IRO) algorithm. 
The IRO algorithm is performed using the subsample of linked reads obtained during the first step of the MTG-Link pipeline.
"""

from __future__ import print_function
import collections
import itertools
import os
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter
from main import gfa_name, subsamplingDir, chunk_size, ext_size
from helpers import Graph
from ProgDynOptim import DynamicMatrixOptim


# Increase the maximum recursion depth in Python.
sys.setrecursionlimit(50000)


#----------------------------------------------------
# index_read function
#----------------------------------------------------
def index_read(read, i, read_rc, seed_size, seedDict):
    """
    To index a read by its seed.

    Args:
        - read: str
            sequence of the current read to index
        - i: str
            position of the current read in readList (list containing all reads' sequences)
        - read_rc: str
            sequence of the reverse complement of the current read
        - seed_size: int
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
        seed = read[:seed_size]
        if seed in seedDict:
            seedDict[seed].append(str(i))
        else:
            seedDict[seed] = [str(i)]

        # Index reverse complement of read by its seed as well.
        seed = read_rc[:seed_size]
        if seed in seedDict:
            seedDict[seed].append("-"+str(i))
        else:
            seedDict[seed] = ["-"+str(i)]

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the function 'index_read()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# find_overlapping_reads function
#----------------------------------------------------
def find_overlapping_reads(assembly, len_read, readList, seed_size, min_overlap, dmax, seedDict):
    """
    To find the reads overlapping with the current assembly's sequence S

    Args:
        - assembly: str
            current assembly's sequence
        - len_read: int
            length of the read from which we want to extend
        - readList: list
            reads of the union
        - seed_size: int
            size of the seed used for indexing the reads
        - min_overlap: int
            minimum overlapping size
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads
        - seedDict: dict
            dictionary 'seedDict' containing reads indexed by their seed
            key = seed's sequence ; value = list of positions of reads having this seed in readList

    Returns:
        - overlapping_reads: list
            list 'overlapping_reads' containing all the overlapping reads' sequences, along with the index of the beginning of the overlap on S, and the index of the beginning of extension on R,
            referenced as [read's sequence, index of beginning of overlap, index of beginning of extension]
            this list is sorted automatically by smallest i, e.g. by largest overlap
    """
    try:
        overlapping_reads = []

        # Get the putative reads (e.g. reads having a seed onto the current assembly's sequence).
        for i in range(len(assembly)-len_read+1, len(assembly)-min_overlap-seed_size):
            seed = assembly[i:i+seed_size]
            if seed in seedDict:
                putative_reads = seedDict[seed]

                # For each putative read, search for an overlap between the current assembly's sequence and the putative read, using the dynamic programmation.
                for put_read in putative_reads:

                    # Get the sequence of the read.
                    if '-' in str(put_read):
                        read = str(Seq(readList[int(put_read.split('-')[1])]).reverse_complement())
                    else:
                        read = readList[int(put_read)]
                    
                    # Perform the alignment with the optimized dynamic programmation.
                    dm = DynamicMatrixOptim(assembly[(i+seed_size):], read[seed_size:], dmax)
                    dist, posR = dm.getEditDistanceAndGenomePosition()
                    
                    # Overlap found.
                    if dist != None:
                        posExt = posR + seed_size
                        overlapping_reads.append([read, i, posExt])

        return overlapping_reads

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the function 'find_overlapping_reads()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# extend function
#----------------------------------------------------
def extend(assembly, len_read, input_seqName, STOP, seedDict, assemblyHash, readList, seed_size, min_overlap, abundance_minList, dmax, max_length, iroLog):
    """
    To extend a read's sequence with overlapping reads.
    This is an iterative function.
    The Boolean value it returns represents the success of the gap-filling
    NB: extGroup is a dictionary containing the extension's sequence as key, and the reads sharing this extension as value
        (value format: [read's sequence, index of beginning of overlap])
    If we use the 'graph' module: def extend(S, read, a, seedDict, graph):

    Args:
        - assembly: str
            current assembly's sequence
        - len_read: int
            length of the read from which we want to extend
        - input_seqName: str
            name of the input k-mers sequences
        - STOP: str
            k-mer STOP (target)
        - seedDict: dict
            dictionary 'seedDict' containing reads indexed by their seed
            key = seed's sequence ; value = list of positions of reads having this seed in readList
        - assemblyHash = hashtable/dict
            hashtable/dictionary indicating if the search for overlapping reads has already been performed on the corresponding sequence (key):
            key = the last 70 bp of the current assembly's sequence ; value = Boolean value (0: overlapping reads search not performed / 1: overlapping reads search performed)
        - readList: list
            reads of the union
        - seed_size: int
            size of the seed used for indexing the reads
        - min_overlap: int
            minimum overlapping size
        - abundance_minList: list
            list of minimal abundance(s) of reads used for gapfilling
            extension's groups having less than this number of reads are discarded from the graph
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads
        - max_length: int
            maximum assembly length (bp)
        - iroLog: file
            the temporary (e.g. incomplete) gap-filled sequences will be saved in this log file

    Return:
        str, Boolean
            - the gap-filled sequence (assembly) and a Boolean variable equal to True if a solution is found (e.g. we arrived to STOP kmer)
            OR
            - the reason why the gap-filling failed and a Boolean variable equal to False if no solution is found
    """
    try:
        # Base cases.
        if STOP in assembly[-len_read:]:
            return assembly, True

        if len(assembly) > max_length:
            return "\n|S| > max_length", False

        if len(assembly) >= 70:
            # Check that we didn't already search for overlapping reads on this region (e.g. on the last 70 bp of the current assembly's sequence).
            if assemblyHash[assembly[-70:]] == 1:
                return "\nPath already explored: No solution", False
                
        # Search for reads overlapping with the current assembly's sequence.
        overlapping_reads = find_overlapping_reads(assembly, len_read, readList, seed_size, min_overlap, dmax, seedDict)
        if not overlapping_reads:
            with open(iroLog, "a") as log:
                log.write("\n>" + input_seqName + " _ No_read_overlapping")
                log.write("\n"+str(assembly)+"\n")
            return "\nNo overlapping reads", False

        # Group the overlapping reads by their extension.
        extGroup = {}

        # Populate extGroup.
        '''NB: overlapping_reads list sorted automatically by smallest i, e.g. by largest overlap'''
        for (read_seq, posS, posExt) in overlapping_reads:

            # If no extension, don't add it to extGroup.
            if read_seq[posExt:] == "":
                continue

            # Add first extension to extGroup.
            if len(extGroup) == 0:
                extGroup[read_seq[posExt:]] = [[read_seq, posS]]

            # Add all extensions to extGroup.
            elif len(extGroup) > 0:
                for extension in extGroup:

                    # Check that current extension is smaller than the one(s) in extGroup.
                    if len(read_seq[posExt:]) < len(extension):
                        # Current extension already in extGroup.
                        if read_seq[posExt:] == extension[:len(read_seq[posExt:])]:
                            new_extension = read_seq[posExt:]
                            extGroup[new_extension] = extGroup[extension]
                            extGroup[new_extension].append([read_seq, posS])
                            del extGroup[extension]
                            added_to_extGroup = True
                            break
                        # Current extension is partially in extGroup.
                        elif read_seq[posExt] == extension[0]:
                            i = 1
                            while i < len(read_seq[posExt:]):
                                if read_seq[posExt+i] == extension[i]:
                                    i += 1
                                else:
                                    break
                            new_extension = extension[:i]
                            extGroup[new_extension] = extGroup[extension]
                            extGroup[new_extension].append([read_seq, posS])
                            del extGroup[extension]
                            added_to_extGroup = True
                            break
                        # Current extension not already in extGroup.
                        else:
                            added_to_extGroup = False

                    # Current extension is not smaller than the one(s) in extGroup.
                    else:
                        # Current extension already in extGroup.
                        if read_seq[posExt:posExt+len(extension)] == extension:
                            extGroup[extension].append([read_seq, posS])
                            added_to_extGroup = True
                            break
                        # Current extension is partially in extGroup.
                        elif read_seq[posExt] == extension[0]:
                            i = 1
                            while i < len(extension):
                                if read_seq[posExt+i] == extension[i]:
                                    i += 1
                                else:
                                    break
                            new_extension = extension[:i]
                            extGroup[new_extension] = extGroup[extension]
                            extGroup[new_extension].append([read_seq, posS])
                            del extGroup[extension]
                            added_to_extGroup = True
                            break
                        # Current extension not already in extGroup.
                        else:
                            added_to_extGroup = False
                            
                # Current extension not already in extGroup.
                if not added_to_extGroup:
                    extGroup[read_seq[posExt:]] = [[read_seq, posS]]

            # Sort extGroup by the smallest extension.
            extGroup = collections.OrderedDict(sorted(extGroup.items(), key=lambda t: len(t[0])))

        # Update 'assemblyHash' to indicate that we performed the search for overlapping reads on this region (e.g. on the last 70 bp of the current assembly's sequence).
        assemblyHash[assembly[-70:]] = 1

        # Filter extGroup by the number of reads sharing an extension (argument 'abundance_min').
        for abundance_min in abundance_minList:
            extGroup_filtered = extGroup.copy()
            for extension in list(extGroup_filtered.keys()):
                if len(extGroup_filtered[extension]) < abundance_min:
                    del extGroup_filtered[extension]
            
            # Iterate over the values of abundance_min only if number of reads sharing an extension < 'abundance_min'.
            if not extGroup_filtered:
                continue
            else:
                break

        # If number of reads sharing an extension < minimal 'abundance_min' provided, stop the extension.
        if not extGroup_filtered:
            with open(iroLog, "a") as log:
                log.write("\n>" + input_seqName + " _ No_extGroup")
                log.write("\n"+str(assembly)+"\n")
            return "\nNo extension", False

        # Sort extGroup by the extension whose read has the largest overlap with the current assembly's sequence (smallest i). 
        '''NB: values of extGroup sorted by reads having the larger overlap'''
        extGroup_filtered = collections.OrderedDict(sorted(extGroup_filtered.items(), key=lambda  t: t[1][0][1]))

        # If there are many possible extensions, iterate over the 5 first extensions.
        if len(extGroup_filtered) > 5:
            final_extGroup = dict(itertools.islice(extGroup_filtered.items(), 5))
        else:
            final_extGroup = extGroup_filtered.copy()

        # Iterative extension of the assembly's sequence S.
        for extension in final_extGroup:

            # Update 'assemblyHash' with the new region for which we will search for overlapping reads (with value '0' if search not already performed, or with value '1' if search already performed).
            if (assembly+extension)[-70:] in assemblyHash.keys():
                assemblyHash[(assembly+extension)[-70:]] = 1
            else:
                assemblyHash[(assembly+extension)[-70:]] = 0
            
            # Extend the updated assembly sequence (e.g. assembly+extension sequence).
            res, success = extend(assembly+extension, len(extGroup_filtered[extension][0][0]), input_seqName, STOP, seedDict, assemblyHash, readList, seed_size, min_overlap, abundance_minList, dmax, max_length, iroLog)

            # If we find a complete gap-filled sequence (e.g. we reach the kmer STOP), return the assembly sequence along with True.
            if success:
                return res, True

        # If we don't find a complete gap-filled sequence, return the reason why the gap-filling was not successful along with False.
        return res, False

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the function 'extend()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# iro_fill function
#----------------------------------------------------
def iro_fill(gap_label, readList, fasta_file, seed_size, min_overlap, abundance_minList, dmax, max_length):
    """
    To execute the IRO algorithm, that is based on on-the-fly computations of read overlaps and iterative extensions of the current assembly sequence.
    The local assembly step is performed between the sequences surrounding the extended gap (e.g. the kmers of the fasta file), in the forward orientation.

    Args:
        - gap_label: str
            label of the gap
        - readList: list
            reads of the union
        - fasta_file: file
            START and STOP k-mers sequences (source and target)
        - seed_size: int
            size of the seed used for indexing the reads
        - min_overlap: int
            minimum overlapping size
        - abundance_minList: list
            list of minimal abundance(s) of reads used for gapfilling
            extension's groups having less than this number of reads are discarded from the graph
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads
        - max_length: int
            maximum assembly length (bp)

    Return:
        - res: str
            - the assembled sequence from the read containing the k-mer START to the read containing the k-mer STOP (assembly) and a Boolean value equal to True
              if a solution is found (e.g. we arrived to STOP kmer)
            OR
            - the reason why the gap-filling failed and a Boolean value equal to False
              if no solution is found
    """
    try:
        # Initiate the four main variables.
        seedDict = {}
        read_with_startList = []
        pos_read_in_readList = 0
        assemblyHash = {}

        # Initiate the log file.
        iroLog = str(gap_label) + "_iro.log"

        # Get the k-mers gap flanking sequences (kmers START and STOP).
        with open(fasta_file, "r") as fasta_input:
            for record in SeqIO.parse(fasta_input, "fasta"):
                if record.id == "start" or "left" in record.description:
                    START = str(record.seq)
                    input_seqName = record.id
                if record.id == "stop" or "right" in record.description:
                    STOP = str(record.seq)
                    if record.id != input_seqName:
                        input_seqName += "-" + record.id

        # Iterate over the reads of 'readList' to obtain the 'seedDict' dictionary and the 'readWithStart' list.
        for read in readList:

            # Get the reverse complement of the current read.
            read_rc = str(Seq(read).reverse_complement())
            # Seed the read and update the 'seedDict' dictionary.
            index_read(read, pos_read_in_readList, read_rc, seed_size, seedDict)

            # Search if the read contains the whole kmer START's sequence (e.g. source k-mer) and update the 'read_with_startList' list.
            if START in read:
                read_with_startList.append([str(pos_read_in_readList), read.index(START)])
            elif START in read_rc:
                read_with_startList.append(["-"+str(pos_read_in_readList), read_rc.index(START)])

            # Increment the position of the current read in 'readList'
            pos_read_in_readList += 1

        # Sort the 'read_with_startList' list by the minimum extension size (e.g. by the maximum index).
        read_with_startList = sorted(read_with_startList, key=itemgetter(1), reverse=True)
        
        # If there is no read containing the kmer START, raise an exception.
        if not read_with_startList:
            print("\n{}: No read in the dataset provided contains the kmer START... Hence, tentative of gapfilling aborted...".format(gap_label))
            sys.exit(1)

        # Extend the reads containing the whole kmer START's sequence.
        for (pos_read, index) in read_with_startList:

            # Get the sequence of the read.
            if '-' in str(pos_read):
                read = str(Seq(readList[int(pos_read.split('-')[1])]).reverse_complement())
            else:
                read = readList[int(pos_read)]

            # Extend the assembly sequence (e.g. the current read containing the whole kmer START's sequence).
            assemblyHash[read[-70:]] = 0
            res, success = extend(read, len(read), input_seqName, STOP, seedDict, assemblyHash, readList, seed_size, min_overlap, abundance_minList, dmax, max_length, iroLog)

            # Case of successful gap-filling.
            break

        return res, success

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the function 'iro_fill()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# iro_assembly function
#----------------------------------------------------
def iro_assembly(gap_label, gap, left_scaffold, right_scaffold, seq_L, seq_R, max_length, seed_size, min_overlap, abundance_minList, dmax):
    """
    To perform the Local Assembly step using an Iterative Read Overlap (IRO) algorithm.  
    The IRO algorithm is based on on-the-fly computations of read overlaps and iterative extensions of the current assembly sequence. This module is executed on the reads of the union.
    This consists of four main steps: get the input FASTA file, get the 'readList' list, perform the local assembly with the IRO algorithm, process the output of the IRO algorithm.

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
        - seed_size: int
            size of the seed used for indexing the reads
        - min_overlap: int
            minimum overlapping size
        - abundance_minList: list
            list of minimal abundance(s) of reads used for gapfilling
            extension's groups having less than this number of reads are discarded from the graph
        - dmax: int
            maximum number of gaps/substitutions allowed in the inexact overlap between reads

    Return:
        - gapfillingFile: file
            file containing the obtained gap-filled sequence(s)
    """
    #----------------------------------------------------
    # Input FASTA file, with offset of size ext removed
    #----------------------------------------------------
    try:
        # Get a FASTA file containing the input sequences (start and stop kmers).
        fastaFile = "{}.{}.g{}.c{}.start_stop.fa".format(gfa_name, str(gap_label), gap.length, chunk_size)
        with open(fastaFile, "w") as fasta:

            # Start sequence (left kmer).
            line1 = ">ctg{}_start _ len_31_bp (left)\n".format(left_scaffold.scaffold)
            line2 = seq_L[(left_scaffold.slen - ext_size - 31):(left_scaffold.slen - ext_size)]
            
            # Stop sequence (right kmer).
            line3 = "\n>ctg{}_stop _ len_31_bp (right)\n".format(right_scaffold.scaffold)
            line4 = seq_R[ext_size:(ext_size + 31)]

            fasta.writelines([line1, line2, line3, line4])

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the input Fasta file (containing START/STOP k-mers) creation")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Save sequences of the read subsampling in a list
    #----------------------------------------------------
    try:
        # Create the list 'readList' containing all reads' sequences extracted from the Read Subsampling step.
        readList = []
        union_readsFile = "{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, str(gap_label), gap.length, chunk_size)
        subreadsFile = os.path.join(subsamplingDir, union_readsFile)
        with open(subreadsFile, "r") as reads_file:
            readList = [str(read.seq) for read in SeqIO.parse(subreadsFile, "fastq")]

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the 'readList' list creation")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # IRO algorithm
    #----------------------------------------------------
    try:
        # Convert the list of abundance min to a string.
        abundance_min = '-'.join(map(str, abundance_minList))
        
        print("\nGap-filling of {} for s={}, o={}, a={} and dmax={} (union)".format(str(gap_label), seed_size, min_overlap, abundance_min, dmax))

        # Determine the maximum assembly length (bp) if the gap length is known.
        ## Add twice the extension size to the gap length (extension on both sides of the gap) and twice the length of reads (e.g. 2x 150 bp) to be large
        if gap.length >= max_length:
            max_length = gap.length + 2*ext_size + 2*150

        # Perform the local assembly step with the IRO algorithm.
        res, success = iro_fill(gap_label, readList, fastaFile, seed_size, min_overlap, abundance_minList, dmax, max_length)

    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the initialization of the IRO algorithm")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Output from IRO algorithm
    #----------------------------------------------------
    try:
        # Output file containing the gap-filled sequence(s)
        gapfillingFile = "{}.{}.g{}.c{}.s{}.o{}.a{}.dmax{}.bxu.insertions.fasta".format(gfa_name, str(gap_label), gap.length, chunk_size, seed_size, min_overlap, abundance_min, dmax)

        # Case of unsuccessful gap-filling.
        if not success:
            print("\n{}: {}".format(gap_label, res))

            # Save the reason why there is no complete gap-filling in a log file.
            iroLog = str(gap_label) + "_iro.log"
            with open(iroLog, "a") as log:
                log.write("\n" + res)

        # Case of successful gap-filling.
        if success:
            print("\n{}: Successful Gap-filling !". format(gap_label))

            # Get the k-mers gap flanking sequences (kmers START and STOP).
            START = seq_L[(left_scaffold.slen - ext_size - 31):(left_scaffold.slen - ext_size)]
            STOP = seq_R[ext_size:(ext_size + 31)]
            input_seqName = "ctg{}_start-ctg{}_stop".format(gap.left, gap.right)

            # Save and pre-process the gap-filled sequence obtained for further qualitative evaluation.
            with open(gapfillingFile, "a") as gapfilling_file:
                assembly_startBegin = res.index(START)
                assembly_stopBegin = res.index(STOP)
                seq = res[assembly_startBegin+len(START):assembly_stopBegin]
                seq_name = "assembly." + input_seqName + " len_" + str(len(seq))
                gapfilling_file.write(">" + seq_name)
                gapfilling_file.write("\n" + seq + "\n")
    
    except Exception as e:
        print("\nFile 'IRO.py': Something wrong with the processing of the output from the IRO algorithm")
        print("Exception-")
        print(e)
        sys.exit(1)


    return gapfillingFile
    
    