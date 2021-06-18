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

"""Module 'ReadSubsampling.py': Read Subsampling

The module 'ReadSubsampling.py' enables to perform the read subsampling step of the MTG-Link gap-filling pipeline.
Linked reads whose barcode is observed in chunk regions surrounding the gap are extracted, and constitute the read subsample used in the local assembly step.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
from main import gfa_name, bamFile, readsFile, indexFile, rbxuFile, chunk_size, barcodes_min_freq


#----------------------------------------------------
# extract_barcodes function
#----------------------------------------------------
def extract_barcodes(bam, gap_label, region, barcodes_occ):
    """
    To extract the barcodes of reads mapping on chunk regions, with `LRez extract`. 
    `LRez extract` enables to extract the list of barcodes from a given region of a BAM file.

    Args:
        - bam: file
            indexed BAM file obtained after mapping the linked reads onto the draft assembly
        - gap_label: str
            label of the gap
        - region: str
            chunk region on which to extract the barcodes
        - barcodes_occ: dict
            this function will output a dictionary 'barcodes_occ' containing the occurences of each barcode
            key = barcode sequence ; value = number of occurences

    Return/Output:
        - barcodes_occ: dict
            dictionary 'barcodes_occ' containing the occurences for each barcode extracted on the chunk region
            key = barcode sequence ; value = number of occurences
    """
    try:
        # LRez extract. 
        ## Parameter -d to include duplicate barcodes
        command = ["LRez", "extract", "--bam", bam, "--region", region, "-d"]
        extractLog = str(gap_label) + "_extract.log"
        tmp_barcodesFile = str(gap_label) + "_extract_stdout.txt"

        with open(tmp_barcodesFile, "w+") as f, open(extractLog, "a") as log:
            subprocess.run(command, stdout=f, stderr=log)
            f.seek(0)

            # Save the barcodes and their occurences in the dict 'barcodes_occ'.
            for line in f.readlines():
                ## Remove the '\n' at the end of the sequence (if any)
                barcode_seq = line.split('\n')[0]
                ## Count occurences of each barcode and save them in the dict 'barcodes_occ'
                if barcode_seq in barcodes_occ:
                    barcodes_occ[barcode_seq] += 1
                else:
                    barcodes_occ[barcode_seq] = 1

        # Remove the raw files obtained from `LRez extract`. 
        subprocess.run(["rm", tmp_barcodesFile])
        if os.path.getsize(extractLog) == 0:
            subprocess.run(["rm", extractLog])

        return barcodes_occ

    except Exception as e:
        print("\nFile 'ReadSubsampling.py': Something wrong with the function 'extract_barcodes()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# get_reads function
#----------------------------------------------------
def get_reads(reads, index, gap_label, barcodes, out_reads):
    """
    To extract the reads associated to the barcodes extracted on chunk regions, with `LRez query fastq`. 
    `LRez query fastq` enables to query a barcodes index and a fastq file to retrieve alignments containing the query barcodes.

    Args:
        - reads: file
            barcoded FASTQ file of linked reads
        - index: file
            index of barcodes
        - gap_label: str
            label of the gap
        - barcodes: file
            barcodes of the union
        - out_reads: file
            this function will output a file containing the reads of the union

    Return/Output:
        - out_reads: dict
            file containing the reads of the union
    """
    try:
        # LRez query fastq.
        command = ["LRez", "query", "fastq", "--fastq", reads, "--index", index, "--list", barcodes]
        getreadsLog = str(gap_label) + "_getreads.log"

        with open(getreadsLog, "a") as log:
            subprocess.run(command, stdout=out_reads, stderr=log)

        # Remove the raw files obtained from `LRez query fastq`.
        if os.path.getsize(getreadsLog) == 0:
            subprocess.run(["rm", getreadsLog])

        return out_reads

    except Exception as e:
        print("\nFile 'ReadSubsampling.py': Something wrong with the function 'get_reads()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# read_subsampling function
#----------------------------------------------------
def read_subsampling(gap_label, gap, left_scaffold, right_scaffold, chunk_L, chunk_R):
    """
    To perform the Read Subsampling step. 
    This step uses the barcode information of the linked read dataset to get a subsample of reads of potential interest for gap-filling. 
    This consists of two main steps: Extract Barcodes, Get Reads.

    Args:
        - gap_label: str
            label of the gap
        - gap: object
            Gap object for the current gap
        - left_scaffold: object
            Scaffold object for the left scaffold
        - right_scaffold: object
            Scaffold object for the right scaffold 
        - chunk_L: int
            size of the left chunk
        - chunk_R: int
            size of the right chunk
    
    Return:
        - union_summary: list
            list 'union_summary' containing the gap ID, the names of the left and right flanking sequences, the gap size, the chunk size, and the number of barcodes and reads extracted on the chunk regions 
            (these reads will further be used to perform the gap-filling)
    """
    #----------------------------------------------------
    # Extract Barcodes
    #----------------------------------------------------    
    try:
        # Initiate a dictionary to count the occurences of each barcode extracted on the chunk regions
        barcodes_occ = {}
        
        # Obtain the left barcodes extracted on the left region (left chunk) and store the barcodes and their occurences in the dict 'barcodes_occ'.
        left_region = left_scaffold.chunk(chunk_L)
        extract_barcodes(bamFile, gap_label, left_region, barcodes_occ)

        # Obtain the right barcodes extracted on the right region (right chunk) and store the barcodes and their occurences in the dict 'barcodes_occ'.
        right_region = right_scaffold.chunk(chunk_R)
        extract_barcodes(bamFile, gap_label, right_region, barcodes_occ)

        # Do the union of the barcodes on both left and right regions (e.g. both left and right chunks). 
        union_barcodesFile = "{}.{}.g{}.c{}.bxu".format(gfa_name, str(gap_label), gap.length, chunk_size)
        with open(union_barcodesFile, "w") as union_barcodes:
            ## Filter barcodes by the minimal frequence of barcodes observed in the union set from the two flanking gap sequences
            for (barcode, occurences) in barcodes_occ.items():
                if occurences >= barcodes_min_freq:
                    union_barcodes.write(barcode + "\n")

    except Exception as e:
        print("\nFile 'ReadSubsampling.py': Something wrong with the extraction of barcodes")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Get Reads
    #----------------------------------------------------
    try:
        # Union: extract the reads associated with the barcodes extracted.
        if rbxuFile == "":
            union_readsFile = "{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, str(gap_label), gap.length, chunk_size)
            with open(union_readsFile, "w") as union_reads:
                get_reads(readsFile, indexFile, gap_label, union_barcodesFile, union_reads)
        
        # If the reads of the union are already extracted, use the corresponding file.
        else:
            union_readsFile = os.path.abspath(rbxuFile)

    except Exception as e:
        print("\nFile 'ReadSubsampling.py': Something wrong with the extraction of reads")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Summary of union (barcodes and reads)
    #----------------------------------------------------
    try:
        bxu = sum(1 for line in open(union_barcodesFile, "r"))
        rbxu = sum(1 for line in open(union_readsFile, "r"))/4
        union_summary = [str(gap.identity), str(gap.left), str(gap.right), gap.length, chunk_size, bxu, rbxu]
        
    except Exception as e:
        print("\nFile 'ReadSubsampling.py': Something wrong with the summary")
        print("Exception-")
        print(e)
        sys.exit(1)


    return union_summary

