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

"""Module 'qualitativeEvaluation.py': Qualitative Evaluation

The module 'qualitativeEvaluation.py' enables to perform a qualitative evaluation of the gap-filled sequence(s) obtained during the Local Assembly step.
NUCmer alignments are performed between the gap-filled sequence(s) and a Reference, and a quality score is calculated based on this alignment.
"""

from __future__ import print_function
import csv
import os
import re
import subprocess
import sys
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO
import main
from helpers import Gap, Scaffold


#----------------------------------------------------
# stats_align function
#----------------------------------------------------
def statsNucmerAlignments(gapLabel, queryFile, referenceFile, extSize, prefix, outputDir):
    """
    To perform statistics on the Nucmer alignment(s) obtained between Reference sequence(s) and Query sequence(s).
    The query sequence(s) are the gap-filled sequence(s).
    The Reference can be either the reference sequence or the gap flanking contigs' sequences.

    Args:
        - gapLabel: str
            label of the gap
        - queryFile: file
            file containing the gap-filled sequence(s) (Query)
        - referenceFile: file
            file containing the reference sequence(s) (Reference)
        - extSize: int
            size of the gap extension on both sides (bp); determine start/end of the local assembly
        - prefix: str
            prefix of the output files
        - outputDir: str
            name of the output directory

    Return/Output:
        files containing the statistics on the Nucmer alignments "Query vs Ref"
    """
    try:
        # Call to script 'stats_alignment.py'.
        scriptPath = sys.path[0]
        stats_alignment_command = os.path.join(scriptPath, "stats_alignment.py")
        command = [stats_alignment_command, "-qry", str(queryFile), "-ref", str(referenceFile), "-ext", str(extSize), "-p", prefix, "-out", outputDir]
        statsLog = str(gapLabel) + "_statsAlignments.log"

        try:
            with open(statsLog, "a") as log:
                subprocess.run(command, stderr=log)
        except IOError as err:
            print("File 'qualitativeEvaluation.py', function 'statsNucmerAlignments()': Unable to open or write to the file {}. \nIOError-{}".format(statsLog, err))
            sys.exit(1)

        # Remove the raw file obtained from `stats_alignment.py`.
        if os.path.getsize(statsLog) <= 0:
            subprocess.run(["rm", statsLog])

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the function 'statsNucmerAlignments()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# getPositionForEdgesOfGFA function
#----------------------------------------------------
def getPositionForEdgesOfGFA(orient1, orient2, length1, length2, extSize):
    """
    To get the beginning and ending positions of the overlap of the current alignment.

    Args:
        - orient1: str
            orientation sign of first sequence
        - orient2: str
            orientation sign of second sequence
        - length1: int
            length of first sequence
        - length2: int
            length of second sequence
        - extSize: int
            size of the gap extension on both sides (bp); determine start/end of the local assembly

    Return:
        - positions: list
            list containing the beginning and ending positions of the overlap
    """
    try:
        # Same orientations.
        if orient1 == orient2:

            # Forward orientations.
            if orient1 == "+":
                beg1 = str(length1 - extSize)
                end1 = str(length1) + "$"  
                beg2 = str(0)
                end2 = str(extSize)

            # Reverse orientations.
            elif orient1 == "-":
                beg1 = str(0)
                end1 = str(extSize)
                beg2 = str(length2 - extSize)
                end2 = str(length2) + "$"

        # Opposite orientations.
        elif orient1 != orient2:

            # First sequence in fwd orientation and second sequence in rev orientation.
            if orient1 == "+":
                beg1 = str(length1 - extSize)
                end1 = str(length1) + "$"
                beg2 = str(length2 - extSize)
                end2 = str(length2) + "$"

            # First sequence in rev orientation and first sequence in fwd orientation.
            elif orient1 == "-":
                beg1 = str(0)
                end1 = str(extSize)
                beg2 = str(0)
                end2 = str(extSize)

        # Get the list 'positions'
        positions = [beg1, end1, beg2, end2]

        if not positions:
            print("File 'qualitativeEvaluation.py, function 'getPositionForEdgesOfGFA()': Unable to create the 'positions' list.", file=sys.stderr)
            sys.exit(1)

        return positions

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the function 'getPositionForEdges()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# getOutputForGfa function
#----------------------------------------------------
def getOutputForGfa(record, extSize, leftName, rightName, leftScaffold, rightScaffold, module):
    """
    To get the ouput variables for updating the GFA when a solution is found for a gap.

    Args:
        - record: object
            Record object of one sequence from the file containing the gap-filled sequence(s)
        - extSize: int
            size of the gap extension on both sides (bp); determine start/end of the local assembly
        - leftName: str
            name of the left scaffold
        - rightName: str
            name of the right scaffold
        - leftScaffold: object
            Scaffold object for the left scaffold
        - rightScaffold: object
            Scaffold object for the right scaffold
        - module: str
            name of the module used for the local assembly step (DBG or IRO)

    Return:
        - outputGFA: list
            list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the gap-filled sequence
    """
    try:
        # Get the information about the sequence.
        sequence = record.seq
        seqLength = len(sequence)
        orientationSign = "+"
        orientation = "fwd"
        quality = record.description.split('Quality ')[1]

        if "bkpt2" in str(record.id):
            orientationSign = "-"
            orientation = "rev"

        # Get the solution name.
        ## Module DBG.
        if module == "DBG":
            solutionNumber = str(record.id).split('_sol_')[1]
            solutionName = str(leftName) +":"+ str(rightName) + "_" + orientation + "_" + solutionNumber
        ## Module IRO.
        if module == "IRO":
            solutionName = str(leftName) +":"+ str(rightName) + "_" + orientation
        solution = solutionName + orientationSign

        if not solution:
            print("File 'qualitativeEvaluation.py, function 'getOutputForGfa()': Unable to create the 'solution' variable.", file=sys.stderr)
            sys.exit(1)

        # Get the position for edges.
        pos_1 = getPositionForEdgesOfGFA(leftScaffold.orient, orientationSign, leftScaffold.slen, seqLength, extSize)
        pos_2 = getPositionForEdgesOfGFA(orientationSign, rightScaffold.orient, seqLength, rightScaffold.slen, extSize)

        # Get the list 'output_for_gfa'.
        outputGFA = [solutionName, seqLength, str(sequence), solution, pos_1, pos_2, quality]

        if not outputGFA:
            print("File 'qualitativeEvaluation.py, function 'getOutputForGfa()': Unable to create the 'outputGFA' list.", file=sys.stderr)
            sys.exit(1)

        return outputGFA

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the function 'getOutputForGfa()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# qualitativeEvaluationOfTheAssembly function
#----------------------------------------------------
def qualitativeEvaluationOfTheAssembly(current_gap, gfaFile, extSize, gapfillingFile, module):
    """
    To perform the Qualitative Evaluation step. 
    This step compares the gap-filled sequence(s) obtained during the Local Assembly step to the Reference.
    The Reference can be either the reference sequence or the gap flanking contigs' sequences.
    This consists of four main steps: Pre-processing of the current gap, getting the Reference File, obtaining the Nucmer Statistics of the alignments, Quality Estimation of gap-filled sequence(s).

    Args:
        - current_gap: str
            current gap identification
        - gfaFile: file
            GFA file containing the gaps' coordinates
        - extSize: int
            size of the gap extension on both sides (bp); determine start/end of the local assembly
        - gapfillingFile: file
            file containing the gap-filled sequence(s)
        - module: str
            name of the module used for the local assembly step (DBG or IRO)
       
    Return:
        - outputGFAList: list of lists
            list of lists, each list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the gap-filled sequence
    """
    #----------------------------------------------------
    # Pre-Processing
    #----------------------------------------------------
    try:
        # Go in the 'outDir' directory.
        try:
            os.chdir(main.outDir)
        except OSError as err:
            print("File 'qualitativeEvaluation.py': Something wrong with specified directory 'outDir'. \nOSError-{}".format(err))
            sys.exit(1)

        # Open the input GFA file.
        gfa = gfapy.Gfa.from_file(gfaFile)
        if not gfa:
            print("File 'qualitativeEvaluation.py', function 'qualitativeEvaluationOfTheAssembly()': Unable to open the input GFA file {}.".format(str(gfaFile)), file=sys.stderr)
            sys.exit(1)
        
        # Get the corresponding Gap line ('G' line).
        for _gap_ in gfa.gaps:
            if str(_gap_) == current_gap:
                current_gap = _gap_
                ## Create the object 'gap' from the class 'Gap'
                gap = Gap(current_gap)
                if not gap:
                    print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to create the object 'gap' from the class 'Gap'.", file=sys.stderr)
                    sys.exit(1)

        # Get some information on the current gap we are working on.
        gapLabel = gap.label()

        # Create two objects ('leftScaffold' and 'rightScaffold') from the class 'Scaffold'.
        leftScaffold = Scaffold(current_gap, gap.left, gfaFile)
        if not leftScaffold:
            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to create the object 'leftScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)
        rightScaffold = Scaffold(current_gap, gap.right, gfaFile)
        if not rightScaffold:
            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to create the object 'rightScaffold' from the class 'Scaffold'.", file=sys.stderr)
            sys.exit(1)

        # Get the gap flanking sequences (e.g. the flanking contigs sequences).
        leftFlankingSeq = str(leftScaffold.sequence())
        if not leftFlankingSeq:
            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to get the left flanking sequence.", file=sys.stderr)
            sys.exit(1)
        rightFlankingSeq = str(rightScaffold.sequence())
        if not rightFlankingSeq:
            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to get the right flanking sequence.", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the 'Pre-Processing' step of the function 'qualitativeEvaluationOfTheAssembly()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #---------------------------------------------------------
    # Reference File (sequence(s) to compare the assembly to)
    #--------------------------------------------------------- 
    try:
        # Qualitative evaluation with the reference sequence ('refDir' provided by the user).
        if main.refDir != "":
            for current_file in os.listdir(main.refDir):
                if str(gapLabel) in current_file:
                    referenceFile = main.refDir +"/"+ str(current_file)
            if not os.path.isfile(referenceFile):
                print("\nWarning: No reference file was found for this gap. The Qualitative Evaluation step will be performed with the gap flanking contigs' sequences.")
                referenceFile = ""

        # Qualitative evalution with the gap flanking contigs' sequences.
        elif (main.refDir == "") or (referenceFile == ""):

            # Merge both left and right flanking contigs sequences into a unique file (referenceFile).
            referenceFile = main.contigDir +"/"+ str(gapLabel) +".g"+ str(gap.length) +".ext"+ str(extSize) + ".contigs.fasta"
            try:
                with open(referenceFile, "w") as ref_fasta:

                    # Left scaffold oriented '+'.
                    if leftScaffold.orient == "+":
                        ref_fasta.write(">" + leftScaffold.name + "_region:" + str(leftScaffold.slen - extSize) + "-" + str(leftScaffold.slen) + "\n")
                        ref_fasta.write(leftFlankingSeq[(leftScaffold.slen - extSize):leftScaffold.slen])
                    # Left scaffold oriented '-' ~ Right scaffold oriented '+'.
                    elif leftScaffold.orient == "-":
                        ref_fasta.write(">" + leftScaffold.name + "_region:0-" + str(extSize) + "\n")
                        ref_fasta.write(str(rc(leftFlankingSeq)[0:extSize]))

                    # Right scaffold oriented '+'.
                    if rightScaffold.orient == "+":
                        ref_fasta.write("\n>" + rightScaffold.name + "_region:0-" + str(extSize) + "\n")
                        ref_fasta.write(rightFlankingSeq[0:extSize])
                    # Right scaffold oriented '-' ~ Left scaffold oriented '+'.
                    elif rightScaffold.orient == "-":
                        ref_fasta.write("\n>" + rightScaffold.name + "_region:" + str(rightScaffold.slen - extSize) + "-" + str(rightScaffold.slen) + "\n")
                        ref_fasta.write(str(rc(rightFlankingSeq)[(rightScaffold.slen - extSize):rightScaffold.slen]))

            except IOError as err:
                print("File 'qualitativeEvaluation.py', function 'qualitativeEvaluationOfTheAssembly()': Unable to open or write to the file {}. \nIOError-{}".format(referenceFile, err))
                sys.exit(1)

        if not os.path.isfile(referenceFile):
            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': The reference file {} doesn't exist.".format(referenceFile), file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the 'Reference File' creation step of the function 'qualitativeEvaluationOfTheAssembly()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Nucmer Statistics: Query vs Reference
    #----------------------------------------------------       
    try:
        # Perform statistics on the Nucmer alignments obtained between the gap-filled sequence(s) (Query) and the reference sequence(s) (Reference). 
        if "/" in gapfillingFile:
            prefix = gapfillingFile.split('/')[-1].split('.bxu')[0]
        else:
            prefix = gapfillingFile.split('.bxu')[0]
        
        statsNucmerAlignments(gapLabel, gapfillingFile, referenceFile, extSize, prefix, main.evalDir)

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the 'Nucmer Statistics' step of the function 'qualitativeEvaluationOfTheAssembly()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    #----------------------------------------------------
    # Quality Estimatation of gap-filled sequence(s)
    #----------------------------------------------------
    try:
        # Reader for the Nucmer alignments' statistics file.
        refQryFile = main.evalDir + "/" + prefix + ".nucmerAlignments.stats"

        if not os.path.exists(refQryFile):
            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': The file containing the Nucmer alignments' statistics {} doesn't exist.".format(str(refQryFile)), file=sys.stderr)
            sys.exit(1)

        else:
            try:
                refQryOutput = open(refQryFile)
            except IOError as err:
                print("File 'qualitativeEvaluation.py', function 'qualitativeEvaluationOfTheAssembly()': Unable to open or write to the file {}. \nIOError-{}".format(str(refQryFile), err))
                sys.exit(1)

            # Local assembly performed with the DBG algorithm.
            if module == "DBG":
                statsReader = csv.DictReader(refQryOutput, \
                                                fieldnames=("Gap", "Gap.Len", "Chunk", "Kmer.Size", "Min.Abundance.Threshold", "Strand", "Solution", "Qry.Len", "Ref", "Ref.Len", \
                                                            "Start.Ref", "End.Ref", "Start.Qry", "End.Qry", "AlignR.Len", "AlignQ.Len", "%.Id", "%.CovR", "%.CovQ", "Frame.R", "Frame.Q", "Quality"), \
                                                delimiter='\t')

            # Local assembly performed with the IRO algorithm.
            if module == "IRO":
                statsReader = csv.DictReader(refQryOutput, \
                                                fieldnames=("Gap", "Gap.Len", "Chunk", "Seed.Size", "Min.Overlap.Size", "Min.Abundance.Threshold", "dmax", "Qry.Len", "Ref", "Ref.Len", \
                                                            "Start.Ref", "End.Ref", "Start.Qry", "End.Qry", "AlignR.Len", "AlignQ.Len", "%.Id", "%.CovR", "%.CovQ", "Frame.R", "Frame.Q", "Quality"), \
                                                delimiter='\t')

            if not statsReader:
                print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to create the reader 'statsReader'.", file=sys.stderr)
                sys.exit(1)

            # Obtain a quality score for each gap-filled sequence.
            outputGFAList = []
            assemblyWithQualityFile = main.assemblyDir +"/"+ gapfillingFile.split('/')[-1].split('.bxu')[0] + ".bxu.insertions_quality.fasta"
            badSolutionsFile = main.outDir + "/bad_solutions.fasta"

            try:
                with open(gapfillingFile, "r") as query, open(assemblyWithQualityFile, "w") as qualified:
                    for record in SeqIO.parse(query, "fasta"):
                        sequence = record.seq

                        ## Local assembly performed with the DBG algorithm.
                        if module == "DBG":
                            if "bkpt1" in record.id:
                                strand = "fwd"
                            elif "bkpt2" in record.id:
                                strand = "rev"

                        ## Local assembly performed with the IRO algorithm.
                        if module == "IRO":
                            recordLabel = (record.id).split("assembly.ctg")[1].split('_START')[0] +"_"+ (record.id).split("-ctg")[1].split('_STOP')[0]
                            strand = "fwd"

                        if not strand:
                            print("File 'qualitativeEvaluation.py, function 'qualitativeEvaluationOfTheAssembly()': Unable to get the strand of the gap-filled sequence {}.".format(str(record.id)), file=sys.stderr)
                            sys.exit(1)

                        #-------------------------------------------------------------------
                        # Ref = Reference Sequence (ex: reference sequence of simulated gap)
                        #-------------------------------------------------------------------
                        if main.refDir != "":
                            try:
                                # Get the list of quality scores for each alignment of one gap-filled sequence against the reference.
                                qualityScoresList = []
                                for row in statsReader:

                                    ## Local assembly performed with the DBG algorithm.
                                    if module == "DBG":
                                        if (row["Solution"] in record.id) and ((strand == "fwd" and row["Strand"] == "fwd") or (strand == "rev" and row["Strand"] == "rev")):
                                            qualityScoresList.append(row["Quality"])

                                    ## Local assembly performed with the IRO algorithm.
                                    if module == "IRO":
                                        if (row["Gap"] == recordLabel):
                                            qualityScoresList.append(row["Quality"])
                                
                                if qualityScoresList == []:
                                    qualityScoresList.append('D')

                                refQryOutput.seek(0)

                                # Global quality score.
                                gapfilledSeqQuality = min(qualityScoresList)
                                
                                record.description = "Quality " + str(gapfilledSeqQuality)
                                SeqIO.write(record, qualified, "fasta")

                                # Update GFA with only the good solutions (the ones having a good quality score).
                                if (len(sequence) > 2*extSize) and (re.match('^.*Quality [AB]$', record.description)):
                                    outputGFA = getOutputForGfa(record, extSize, gap.left, gap.right, leftScaffold, rightScaffold, module)
                                    outputGFAList.append(outputGFA)

                                # Add the bad solutions to a FASTA file containing all bad solutions.
                                else:
                                    #outputGFAList = []     #Pbm if several solutions output
                                    try:
                                        with open(badSolutionsFile, "a") as badSolFile:
                                            badSolFile.write("\n>" + str(gap.left) +":"+ str(gap.right) + "_gf." + str(strand) + " _ len." + str(len(sequence)) + "_qual." + str(gapfilledSeqQuality) +"\n")
                                            badSolFile.write(str(sequence))
                                    except IOError as err:
                                        print("File 'qualitativeEvaluation.py', function 'qualitativeEvaluationOfTheAssembly()': Unable to open or write to the file {}. \nIOError-{}".format(str(badSolutionsFile), err))
                                        sys.exit(1)
                            
                            except Exception as e:
                                print("File 'qualitativeEvaluation.py': Something wrong with the case 'Ref = Reference sequence' of the 'Quality Estimation' step of the function 'qualitativeEvaluationOfTheAssembly()'")
                                print("Exception-{}".format(e))
                                sys.exit(1)

                        #----------------------------------------------------
                        # Ref = Gap Flanking Contigs' Sequences
                        #----------------------------------------------------
                        else:
                            try:
                                # Get the list of quality scores for each alignment of one gap-filled sequence against the reference.
                                qualityLeftScoresList = []
                                qualityRightScoresList = []
                                for row in statsReader:

                                    ## Local assembly performed with the DBG algorithm.
                                    if module == "DBG":
                                        if (row["Solution"] in record.id) and ((strand == "fwd" and row["Strand"] == "fwd") or (strand == "rev" and row["Strand"] == "rev")) and (row["Ref"] == leftScaffold.name):
                                            qualityLeftScoresList.append(row["Quality"])
                                        elif (row["Solution"] in record.id) and ((strand == "fwd" and row["Strand"] == "fwd") or (strand == "rev" and row["Strand"] == "rev")) and (row["Ref"] == rightScaffold.name):
                                            qualityRightScoresList.append(row["Quality"])

                                    ## Local assembly performed with the IRO algorithm.
                                    if module == "IRO":
                                        if (row["Gap"] == recordLabel) and (row["Ref"] == leftScaffold.name):
                                            qualityLeftScoresList.append(row["Quality"])
                                        elif (row["Gap"] == recordLabel) and (row["Ref"] == rightScaffold.name):
                                            qualityRightScoresList.append(row["Quality"])

                                if qualityLeftScoresList == []:
                                    qualityLeftScoresList.append('D')
                                if qualityRightScoresList == []:
                                    qualityRightScoresList.append('D')

                                refQryOutput.seek(0)

                                # Global quality score.
                                gapfilledSeqQuality = min(qualityLeftScoresList) + min(qualityRightScoresList)

                                record.description = "Quality " + str(gapfilledSeqQuality)
                                SeqIO.write(record, qualified, "fasta")

                                # Update GFA with only the good solutions (the ones having a good quality score).
                                if (len(sequence) > 2*extSize) and (re.match('^.*Quality [AB]{2}$', record.description)):
                                    outputGFA = getOutputForGfa(record, extSize, gap.left, gap.right, leftScaffold, rightScaffold, module)
                                    outputGFAList.append(outputGFA)

                                # Add the bad solutions to a FASTA file containing all bad solutions.
                                else:
                                    #outputGFAList = []         #Pbm if several solutions output
                                    try:
                                        with open(badSolutionsFile, "a") as badSolFile:
                                            badSolFile.write("\n>" + str(gap.left) +":"+ str(gap.right) + "_gf." + str(strand) + " _ len." + str(len(sequence)) + "_qual." + str(gapfilledSeqQuality) +"\n")
                                            badSolFile.write(str(sequence))
                                    except IOError as err:
                                        print("File 'qualitativeEvaluation.py', function 'qualitativeEvaluationOfTheAssembly()': Unable to open or write to the file {}. \nIOError-{}".format(str(badSolutionsFile), err))
                                        sys.exit(1)

                            except Exception as e:
                                print("File 'qualitativeEvaluation.py': Something wrong with the case 'Ref = Gap Flanking Contigs' Sequences' of the 'Quality Estimation' step of the function 'qualitativeEvaluationOfTheAssembly()'")
                                print("Exception-{}".format(e))
                                sys.exit(1)

                    qualified.seek(0)

            except IOError as err:
                print("File 'qualitativeEvaluation.py', function 'qualitativeEvaluationOfTheAssembly()': Unable to open the file {} or to open/write to the file {}. \nIOError-{}".format(str(gapfillingFile), str(assemblyWithQualityFile), err))
                sys.exit(1)

            # Remove the 'gapfillingFile' once done with it.
            subprocess.run(["rm", gapfillingFile])

        # Close the 'refQryOutput' object (which opened the 'refQryFile' file).
        refQryOutput.close()

    except Exception as e:
        print("File 'qualitativeEvaluation.py': Something wrong with the 'Quality Estimation' step of the function 'qualitativeEvaluationOfTheAssembly()'")
        print("Exception-{}".format(e))
        sys.exit(1)


    return outputGFAList

