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

"""Module 'QualEval.py': Qualitative Evaluation

The module 'QualEval.py' enables to perform a qualitative evaluation of the gap-filled sequence(s) obtained during the local assembly step.
NUCmer alignments are performed between the gap-filled sequence(s) and a Reference, and a quality score is calculated based on this alignment.
"""

from __future__ import print_function
import csv
import os
import re
import subprocess
import sys
from gfapy.sequence import rc
from Bio import SeqIO
from main import refDir, outDir, assemblyDir, contigDir, evalDir, ext_size


#----------------------------------------------------
# stats_align function
#----------------------------------------------------
def stats_align(gap_label, qry_file, ref_file, ext, prefix, out_dir):
    """
    To perform statistics on the alignment between a reference sequence and query sequence(s).
    The query sequence(s) are the gap-filled sequence(s).
    The Reference can be either the reference sequence or the gap flanking contigs.

    Args:
        - gap_label: str
            label of the gap
        - qry_file: file
            file containing the gap-filled sequence(s) (Query)
        - ref_file: file
            file containing the reference sequence(s) (Reference)
        - ext: int
            size of the extension on both gap sides
        - prefix: str
            prefix of the output files
        - out_dir: str
            name of the output directory

    Return/Output:
        Files containing the statistics on the alignments Query vs Ref
    """
    try:
        # `stats_alignment.py`.
        scriptPath = sys.path[0]
        stats_align_command = os.path.join(scriptPath, "stats_alignment.py")
        command = [stats_align_command, "-qry", qry_file, "-ref", ref_file, "-ext", ext, "-p", prefix, "-out", out_dir]
        statsLog = str(gap_label) + "_stats_align.log"

        with open(statsLog, "a") as log:
            subprocess.run(command, stderr=log)

        # Remove the raw file obtained from `stats_alignment.py`.
        if os.path.getsize(statsLog) <= 0:
            subprocess.run(["rm", statsLog])

    except Exception as e:
        print("\nFile 'QualEval.py': Something wrong with the function 'stats_align()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# get_position_for_edges function
#----------------------------------------------------
def get_position_for_edges(orient1, orient2, length1, length2, ext):
    """
    To get the beginning and ending positions of the overlap

    Args:
        - orient1: str
            orientation of first sequence
        - orient2: str
            orientation of second sequence
        - length1: int
            length of first sequence
        - length2: int
            length of second sequence
        - ext: int
            size of the extension on both gap sides

    Return:
        - positions: list
            list containing the beginning and ending positions of the overlap
    """
    try:
        # Same orientations.
        if orient1 == orient2:

            # Forward orientations.
            if orient1 == "+":
                beg1 = str(length1 - ext)
                end1 = str(length1) + "$"  
                beg2 = str(0)
                end2 = str(ext)

            # Reverse orientations.
            elif orient1 == "-":
                beg1 = str(0)
                end1 = str(ext)
                beg2 = str(length2 - ext)
                end2 = str(length2) + "$"

        # Opposite orientations.
        elif orient1 != orient2:

            # First seq in fwd orientation and second seq in rev orientation.
            if orient1 == "+":
                beg1 = str(length1 - ext)
                end1 = str(length1) + "$"
                beg2 = str(length2 - ext)
                end2 = str(length2) + "$"

            # First seq in rev orientation and first seq in fwd orientation.
            elif orient1 == "-":
                beg1 = str(0)
                end1 = str(ext)
                beg2 = str(0)
                end2 = str(ext)

        # Get the list 'positions'
        positions = [beg1, end1, beg2, end2]

        return positions

    except Exception as e:
        print("\nFile 'QualEval.py': Something wrong with the function 'get_position_for_edges()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# get_output_for_gfa function
#----------------------------------------------------
def get_output_for_gfa(record, ext, s1, s2, left_scaffold, right_scaffold):
    """
    To get the ouput variables for updating the GFA when a solution is found for a gap

    Args:
        - record: object
            Record object of one sequence from the file containing the gap-filled sequence(s)
        - ext: int
            size of the extension on both gap sides
        - s1: str
            name of the left scaffold
        - s2: str
            name of the right scaffold
        - left_scaffold: object
            Scaffold object for the left scaffold
        - right_scaffold: object
            Scaffold object for the right scaffold

    Return:
        - output_for_gfa: list
            list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the gap-filled sequence
    """
    try:
        # Get the information about the sequence.
        seq = record.seq
        length_seq = len(seq)
        orient_sign = "+"
        orient = "fwd"
        quality = record.description.split('Quality ')[1]

        if "bkpt2" in str(record.id):
            orient_sign = "-"
            orient = "rev"

        # Get the solution name.
        sol_name = str(s1) +":"+ str(s2) + "_gf_"+ orient
        solution = sol_name + orient_sign

        # Get the position for edges.
        pos_1 = get_position_for_edges(left_scaffold.orient, orient_sign, left_scaffold.slen, length_seq, ext)
        pos_2 = get_position_for_edges(orient_sign, right_scaffold.orient, length_seq, right_scaffold.slen, ext)

        # Get the list 'output_for_gfa'.
        output_for_gfa = [sol_name, length_seq, str(seq), solution, pos_1, pos_2, quality]

        return output_for_gfa

    except Exception as e:
        print("\nFile 'QualEval.py': Something wrong with the function 'get_output_for_gfa()'")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# qual_eval function
#----------------------------------------------------
def qual_eval(gap_label, gap, left_scaffold, right_scaffold, seq_L, seq_R, gapfillingFile):
    """
    To perform the Qualitative Evaluation step. 
    This step compares the gap-filled sequence(s) obtained during the local assembly step to the Reference.
    The Reference can be either the reference sequence or the gap flanking contigs.
    This consists of three main steps: get the Reference, compare the gap-filled sequence(s) to the Reference using Nucmer, estimate the quality of gap-filled sequence(s).

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
        - gapfillingFile: file
            file containing the gap-filled sequence(s)
       
    Return:
        - output_for_gfa: list
            list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the sequence
    """
    #----------------------------------------------------
    # Get the Reference
    #----------------------------------------------------   
    try:
        # Qualitative evaluation with the reference sequence.
        if refDir != "":
            for file_ in os.listdir(refDir):
                if str(gap_label) in file_:
                    refFile = refDir +"/"+ str(file_)
            if not os.path.isfile(refFile):
                print("\nWarning: No reference file was found for this gap. The qualitative evaluation will be performed with the flanking contigs information.")
                refFile = ""

        # Qualitative evalution with the flanking contigs information.
        elif (refDir == "") or (refFile == ""):

            # Merge both left and right flanking contigs sequences into a unique file (refFile).
            refFile = contigDir +"/"+ str(gap_label) +".g"+ str(gap.length) + ".contigs.fasta"
            with open(refFile, "w") as ref_fasta:

                # Left scaffold oriented '+'.
                if left_scaffold.orient == "+":
                    ref_fasta.write(">" + left_scaffold.name + "_region:" + str(left_scaffold.slen-ext_size) + "-" + str(left_scaffold.slen) + "\n")
                    ref_fasta.write(seq_L[(left_scaffold.slen - ext_size):left_scaffold.slen])
                # Left scaffold oriented '-' ~ Right scaffold oriented '+'.
                elif left_scaffold.orient == "-":
                    ref_fasta.write(">" + left_scaffold.name + "_region:0-" + str(ext_size) + "\n")
                    ref_fasta.write(str(rc(seq_L)[0:ext_size]))

                # Right scaffold oriented '+'.
                if right_scaffold.orient == "+":
                    ref_fasta.write("\n>" + right_scaffold.name + "_region:0-" + str(ext_size) + "\n")
                    ref_fasta.write(seq_R[0:ext_size])
                # Right scaffold oriented '-' ~ Left scaffold oriented '+'.
                elif right_scaffold.orient == "-":
                    ref_fasta.write("\n>" + right_scaffold.name + "_region:" + str(right_scaffold.slen-ext_size) + "-" + str(right_scaffold.slen) + "\n")
                    ref_fasta.write(str(rc(seq_R)[(right_scaffold.slen - ext_size):right_scaffold.slen]))

        if not os.path.isfile(refFile):
            print("\nWarning: Something wrong with the specified reference file. Exception-", sys.exc_info())
            sys.exit(1)

    except Exception as e:
        print("\nFile 'QualEval.py': Something wrong with the reference file")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Query vs Reference: Nucmer statistics
    #----------------------------------------------------       
    try:
        # Perform statistics on the Nucmer alignments obtained between the gap-filled sequence(s) (Query) and the reference sequence(s) (Reference). 
        if "/" in gapfillingFile:
            prefix = gapfillingFile.split('/')[-1].split('.bxu')[0]
        else:
            prefix = gapfillingFile.split('.bxu')[0]
        
        stats_align(gap_label, gapfillingFile, refFile, ext_size, prefix, evalDir)

    except Exception as e:
        print("\nFile 'QualEval.py': Something wrong with the initialization of the Nucmer alignments")
        print("Exception-")
        print(e)
        sys.exit(1)


    #----------------------------------------------------
    # Estimate quality of gap-filled sequences
    #----------------------------------------------------  
    try:
        # Reader for alignment stats' files.
        ref_qry_file = evalDir + "/" + prefix + ".ref_qry.alignment.stats"

        if not os.path.exists(ref_qry_file):
            print("\nWarning: The '{}' file doesn't exits".format(ref_qry_file))
            sys.exit(1)

        else:
            ref_qry_output = open(ref_qry_file)

            # Local assembly performed with the DBG algorithm.
            if ".k" in gapfillingFile.split('/')[-1]:
                reader_ref_stats = csv.DictReader(ref_qry_output, \
                                                fieldnames=("Gap", "Len_gap", "Chunk", "Kmer_size", "Abundance_min", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                                                            "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
                                                delimiter='\t')

            # Local assembly performed with the IRO algorithm.
            if ".dmax" in gapfillingFile.split('/')[-1]:
                reader_ref_stats = csv.DictReader(ref_qry_output, \
                                                fieldnames=("Gap", "Len_gap", "Chunk", "Seed_size", "Min_overlap", "Abundance_min", "dmax", "Len_Q", "Ref", "Len_R", \
                                                            "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
                                                delimiter='\t')


            # Obtain a quality score for each gapfilled seq.
            output_for_gfa = []
            assembly_quality_file = os.path.abspath(assemblyDir +"/"+ gapfillingFile.split('.bxu')[0] + ".bxu.insertions_quality.fasta")
            bad_solutions_file = os.path.abspath(outDir + "/bad_solutions.fasta")

            with open(gapfillingFile, "r") as query, open(assembly_quality_file, "w") as qualified:
                for record in SeqIO.parse(query, "fasta"):
                    seq = record.seq

                    ## Local assembly performed with the IRO algorithm
                    if ".k" in gapfillingFile.split('/')[-1]:
                        strand = str(record.id).split('_')[0][-1]

                    ## Local assembly performed with the IRO algorithm
                    if ".dmax" in gapfillingFile.split('/')[-1]:
                        record_label = (record.id).split("assembly.ctg")[1].split('_start')[0] +"_"+ (record.id).split("-ctg")[1].split('_stop')[0]
                        strand = "fwd"

                    #----------------------------------------------------
                    #Ref = reference sequence (ex: of simulated gap)
                    #----------------------------------------------------
                    if refDir != "":
                        # Quality score for stats about the alignment Query vs Ref.
                        quality_ref = []
                        for row in reader_ref_stats:

                            ## Local assembly performed with the DBG algorithm
                            if ".k" in gapfillingFile.split('/')[-1]:
                                if (row["Solution"] in record.id) and (("bkpt1" in record.id and row["Strand"] == "fwd") or ("bkpt2" in record.id and row["Strand"] == "rev")):
                                    quality_ref.append(row["Quality"])

                            ## Local assembly performed with the IRO algorithm
                            if ".dmax" in gapfillingFile.split('/')[-1]:
                                if (row["Gap"] == record_label):
                                    quality_ref.append(row["Quality"])
                        
                        if quality_ref == []:
                            quality_ref.append('D')

                        ref_qry_output.seek(0)

                        # Global quality score.
                        quality_gapfilled_seq = min(quality_ref)
                        
                        record.description = "Quality " + str(quality_gapfilled_seq)
                        SeqIO.write(record, qualified, "fasta")

                        # Update GFA with only the good solutions (the ones having a good quality score).
                        if (len(seq) > 2*ext_size) and (re.match('^.*Quality [AB]$', record.description)):
                            gfa_output = get_output_for_gfa(record, ext_size, gap.left, gap.right, left_scaffold, right_scaffold)
                            output_for_gfa.append(gfa_output)

                        # Add the bad solutions to a FASTA file containing all bad solutions.
                        else:
                            #output_for_gfa = []     #Pbm if several solutions output
                            with open(bad_solutions_file, "a") as bad_file:
                                bad_file.write("\n>" + str(gap.left) +":"+ str(gap.right) + "_gf_" + str(strand) + " _ len_" + str(len(seq)) + "_qual_" + str(quality_gapfilled_seq) +"\n")
                                bad_file.write(str(seq))

                    #----------------------------------------------------
                    #Ref = flanking contigs' sequences
                    #----------------------------------------------------
                    else:
                        # Quality score for stats about the alignments Query vs Extensions.
                        quality_ext_left = []
                        quality_ext_right = []
                        for row in reader_ref_stats:

                            ## Local assembly performed with the DBG algorithm
                            if ".k" in gapfillingFile.split('/')[-1]:
                                if (row["Solution"] in record.id) and (("bkpt1" in record.id and row["Strand"] == "fwd") or ("bkpt2" in record.id and row["Strand"] == "rev")) and (row["Ref"] == left_scaffold.name):
                                    quality_ext_left.append(row["Quality"])
                                elif (row["Solution"] in record.id) and (("bkpt1" in record.id and row["Strand"] == "fwd") or ("bkpt2" in record.id and row["Strand"] == "rev")) and (row["Ref"] == right_scaffold.name):
                                    quality_ext_right.append(row["Quality"])

                            ## Local assembly performed with the IRO algorithm
                            if ".dmax" in gapfillingFile.split('/')[-1]:
                                if (row["Gap"] == record_label) and (row["Ref"] == left_scaffold.name):
                                    quality_ext_left.append(row["Quality"])
                                elif (row["Gap"] == record_label) and (row["Ref"] == right_scaffold.name):
                                    quality_ext_right.append(row["Quality"])

                        if quality_ext_left == []:
                            quality_ext_left.append('D')
                        if quality_ext_right == []:
                            quality_ext_right.append('D')

                        ref_qry_output.seek(0)

                        # Global quality score.
                        quality_gapfilled_seq = min(quality_ext_left) + min(quality_ext_right)

                        record.description = "Quality " + str(quality_gapfilled_seq)
                        SeqIO.write(record, qualified, "fasta")

                        #Update GFA with only the good solutions (the ones having a good quality score).
                        if (len(seq) > 2*ext_size) and (re.match('^.*Quality [AB]{2}$', record.description)):
                            gfa_output = get_output_for_gfa(record, ext_size, gap.left, gap.right, left_scaffold, right_scaffold)
                            output_for_gfa.append(gfa_output)

                        #Add the bad solutions to a FASTA file containing all bad solutions.
                        else:
                            #output_for_gfa = []         #Pbm if several solutions output
                            with open(bad_solutions_file, "a") as bad_file:
                                bad_file.write("\n>" + str(gap.left) +":"+ str(gap.right) + "_gf_" + str(strand) + " _ len_" + str(len(seq)) + "_qual_" + str(quality_gapfilled_seq) +"\n")
                                bad_file.write(str(seq))


                qualified.seek(0)

            # If local assembly performed with the DBG algorithm
            if ".k" in gapfillingFile.split('/')[-1]:
                # Remove the 'gapfillingFile' once done with it.
                subprocess.run(["rm", gapfillingFile])

    except Exception as e:
        print("\nFile 'QualEval.py': Something wrong with the estimation of the quality of gap-filled sequence(s)")
        print("Exception-")
        print(e)
        sys.exit(1)


    return output_for_gfa

