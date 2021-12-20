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

"""Module 'mtglink.py': Gap-filling Pipeline MTG-Link

The module 'mtglink.py' enables to process the input and output GFA files, and to process all gaps in a multi-threading way.
It output as well a summary of the gap-filling results.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys
from pathos.multiprocessing import ProcessingPool as Pool
#from multiprocessing import Pool
import gfapy
from Bio import SeqIO
import main
from barcodesExtraction import extractBarcodesFromChunkRegions
from readsRetrieval import retrieveReadsWithLRezQueryFastq
from gapFilling import fillGapByLocalAssemblyAndQualitativeEvaluation
from helpers import updateGFAWithSolution


#----------------------------------------------------
# GFA Pre-Processing
#----------------------------------------------------
try:
    print("\nSTEP 1/3: GFA Pre-Processing")

    # Open the input GFA file.
    gfa = gfapy.Gfa.from_file(main.gfaFile)
    if not gfa:
        print("File 'mtglink.py': Unable to open the input GFA file {}.".format(str(main.gfaFile)), file=sys.stderr)
        sys.exit(1)

    # Go in the 'outDir' directory.
    try:
        os.chdir(main.outDir)
    except OSError as err:
        print("File 'mtglink.py': Something wrong with specified directory 'outDir'. \nOSError-{}".format(err))
        sys.exit(1)
    
    # Create the output GFA file.
    outputGFAFile = main.outDir + "/" + str(main.gfa_name).split('.gfa')[0] + "_mtglink.gfa"

    # Case no gap in the input GFA: if no gap, rewrite all the lines into GFA output.
    if len(gfa.gaps) == 0:
        try: 
            with open(outputGFAFile, "w") as f:
                out_gfa = gfapy.Gfa()
                for line in gfa.lines:
                    out_gfa.add_line(str(line))
                out_gfa.to_file(outputGFAFile)
        except IOError as err:
            print("File 'mtglink.py': Unable to open or write to the output GFA File {}. \nIOError-{}".format(str(outputGFAFile), err))
            sys.exit(1)

    # Case gaps to fill in the input GFA
    gaps = []

    ## If gap and analysis hasn't started, rewrite the H and S lines into GFA output.
    if main.line_gfa == "":
        try:
            with open(outputGFAFile, "w") as f:
                out_gfa = gfapy.Gfa()
                out_gfa.add_line("H\tVN:Z:2.0")
                for line in gfa.segments:
                    out_gfa.add_line(str(line))
                out_gfa.to_file(outputGFAFile)
        except IOError as err:
            print("File 'mtglink.py': Unable to open or write to the output GFA File {}. \nIOError-{}".format(str(outputGFAFile), err))
            sys.exit(1)
        # Convert Gfapy gap line to a string to be able to use it with multiprocessing.
        for _gap_ in gfa.gaps:
            _gap_ = str(_gap_)
            gaps.append(_gap_)
        if len(gaps) == 0:
            print("File 'mtglink.py': Error while converting the Gfapy gap line to a string and appending the 'gaps' list.", file=sys.stderr)
            sys.exit(1)
        
    ## If gap and analysis has started (e.g. '-line' argument provided), start analysis from this line in GFA file input.
    if main.line_gfa != "":
        for _gap_ in gfa.gaps[(main.line_gfa - (len(gfa.segments)+2)):]:
            _gap_ = str(_gap_)
            gaps.append(_gap_)
        if len(gaps) == 0:
            print("File 'mtglink.py': Error while converting the Gfapy gap line to a string and appending the 'gaps' list.", file=sys.stderr)
            sys.exit(1)

except Exception as e:
    print("File 'mtglink.py': Something wrong with the 'GFA Pre-Processing' step.")
    print("Exception-{}".format(e))
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


#----------------------------------------------------   
# Read Subsampling
#----------------------------------------------------
try:
    print("\nSTEP 2/3: Read Subsampling")

    # Go in the 'subsamplingDir' directory.
    try:
        os.chdir(main.subsamplingDir)
    except OSError as err:
        print("File 'mtglink.py': Something wrong with specified directory 'subsamplingDir'. \nOSError-{}".format(err))
        sys.exit(1)

    # For each gap, get the list of barcodes of potential interest (e.g. barcodes of reads mapping on chunk regions) and append the file containing all the lists obtained for all gaps ('allBarcodesFile').
    print("\nSTEP 2a: Barcodes Extraction\n")
    allBarcodesFile = "{}.c{}.f{}.allBarcodesFiles.bx".format(main.gfa_name, main.chunkSize, main.barcodesMinFreq)
    for gap in gaps:
        unionBarcodesFile = extractBarcodesFromChunkRegions(gap, main.gfaFile, main.bamFile, main.chunkSize, main.barcodesMinFreq)
        try:
            with open(allBarcodesFile, "a") as barcFile:
                barcFile.write(str(unionBarcodesFile) + "\n")
        except IOError as err:
            print("File 'mtglink.py': Unable to open or write to the output file containing all the lists of barcodes {}. \nIOError-{}".format(str(allBarcodesFile), err))
            sys.exit(1)

    # For each gap, query the barcodes index and the FASTQ file to retrieve all reads whose barcode is present in the list of barcodes of potential interest.
    print("\nSTEP 2b: Reads Retrieval")
    retrieveReadsWithLRezQueryFastq(main.gfa_name, main.readsFile, main.indexFile, allBarcodesFile, main.threads)

    # For each gap, generate a summary of the number of barcodes and reads extracted from the union of both gap flanking regions.
    print("\nSTEP 2c: Read Subsampling Summary")
    readSubsamplingSummaryFile = "{}.c{}.f{}.readSubsampling_summary.txt".format(main.gfa_name, main.chunkSize, main.barcodesMinFreq)
    try:
        with open(readSubsamplingSummaryFile, "w") as summaryFile:
            legend = ["Left_Scaffold", "Right_Scaffold", "Gap_Size", "Chunk_Size", "MinFreq_Barcodes", "Nb_Barcodes", "Nb_Reads"]
            summaryFile.write('\t'.join(j for j in legend))
            summaryDict = {}

            # Iterate over the gaps ".bxu" and ".bxu.fastq" files to the the summary values.
            for current_file in os.listdir(main.subsamplingDir):

                # Get the number of barcodes.
                if current_file.endswith('.bxu'):
                    bxu = sum(1 for line in open(current_file, "r"))
                    suffixBarcodesFile = str(current_file).split('/')[-1].split('.gfa.')[1]

                    # Get the gap flanking scaffolds name.
                    FlankingScaffold = re.split('\.g.*\.c.*\.f.*\.bxu$', str(suffixBarcodesFile))[0]

                    ## Get the name of the left flanking scaffold.
                    LeftFlankingScaffold_wo_sign = re.split('\+_|\-_', str(suffixBarcodesFile))[0]
                    r1 = re.findall(r"\+_|\-_",str(suffixBarcodesFile))[0]
                    LeftFlankingScaffold_sign = re.split("_", str(r1))[0]
                    LeftFlankingScaffold = LeftFlankingScaffold_wo_sign + LeftFlankingScaffold_sign

                    ## Get the name of the right flanking scaffold.
                    RightFlankingScaffold_suffix = re.split('\+_|\-_', str(suffixBarcodesFile))[1]
                    RightFlankingScaffold_wo_sign = re.split('\+\.|\-\.', str(RightFlankingScaffold_suffix))[0]
                    r2 = re.findall(r"\+\.|\-\.",str(suffixBarcodesFile))[0]
                    RightFlankingScaffold_sign = re.split("\.", str(r2))[0]
                    RightFlankingScaffold = RightFlankingScaffold_wo_sign + RightFlankingScaffold_sign

                    # Get the parameters values ('Gap_Size', 'Chunk_Size', 'MinFreq_Barcodes').
                    suffixParametersValues = re.findall('\.g.*\.c.*\.f.*\.bxu$', str(suffixBarcodesFile))
                    gapSizeValue = suffixParametersValues[0].split('.g')[1].split('.c')[0]
                    chunkSizeValue = suffixParametersValues[0].split('.c')[1].split('.f')[0]
                    minFreqBarcodesValue = suffixParametersValues[0].split('.f')[1].split('.bxu')[0]

                    # Append the 'summaryList' with a list of the corresponding values for one gap.
                    summaryDict[FlankingScaffold] = [LeftFlankingScaffold, RightFlankingScaffold, gapSizeValue, chunkSizeValue, minFreqBarcodesValue, bxu]

                # Get the number of reads.
                elif current_file.endswith('.bxu.fastq'):
                    rbxu = sum(1 for line in open(current_file, "r"))/4
                    suffixReadsFile = str(current_file).split('/')[-1].split('.gfa.')[1]

                    # Get the gap flanking scaffolds name.
                    FlankingScaffold = re.split('\.g.*\.c.*\.f.*\.bxu.fastq$', str(suffixReadsFile))[0]

                    # Append the 'summaryList' with a corresponding reads' value for one gap.
                    if (FlankingScaffold in summaryDict):
                        summaryDict[FlankingScaffold].append(rbxu)

            # Update the 'summaryFile' with the corresponding values for each gap.
            for FlankingScaffold, values in summaryDict.items():
                summaryFile.write("\n" + '\t'.join(str(i) for i in values))

    except IOError as err:
        print("File 'mtglink.py': Unable to open or write to the output summary file of the 'Read Subsampling' step {}. \nIOError-{}".format(str(readSubsamplingSummaryFile), err))
        sys.exit(1)

except Exception as e:
    print("File 'mtglink.py': Something wrong with the 'Read Subsampling' step.")
    print("Exception-{}".format(e))
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


#--------------------------------------------------------------
# Gap-filling (Local Assembly and Qualitative Evaluation steps)
#--------------------------------------------------------------
try:
    print("\nSTEP 3/3: Gap-filling (Local Assembly and Qualitative Evaluation steps)\n")

    # Start multiprocessing.
    p = Pool()

    # Initiate the Boolean variable 'successful_solution' to False.
    gapfillSeqFileExist = False

    for outputGFAList in p.map(fillGapByLocalAssemblyAndQualitativeEvaluation, gaps):
        # Output the 'outputGFAList' results (obtained for each gap) from the function 'fillGapByLocalAssemblyAndQualitativeEvaluation' in the output GFA file.
        print("\nCreating the output GFA file...")

        ## Solution found for the current gap.
        if len(outputGFAList[0]) > 1:          
            for outputGFA in outputGFAList:
                gapfillSeqFile = updateGFAWithSolution(main.outDir, main.gfa_name, outputGFA, outputGFAFile)
                gapfillSeqFileExist = True

        ## No solution found for the current gap.
        else:                                   
            out_gfa = gfapy.Gfa.from_file(outputGFAFile)
            out_gfa.add_line(outputGFAList[0][0])
            out_gfa.to_file(outputGFAFile)

    p.close()

except Exception as e:
    print("File 'mtglink.py': Something wrong with the 'Gap-filling' step.")
    print("Exception-{}".format(e))
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


#----------------------------------------------------
# Files Locations
#----------------------------------------------------
print("\n\nThe results from MTG-Link are saved in: " + main.outDir)
print("The results from the 'Read Subsampling' step are saved in " + main.subsamplingDir)
print("The results from the 'Local Assembly' step are saved in " + main.assemblyDir)
print("The results from the 'Qualitative Evaluation' step are saved in " + main.evalDir)

print("\nSummary of the union: " + main.subsamplingDir +"/"+ readSubsamplingSummaryFile)
print("GFA output file: " + outputGFAFile)
if gapfillSeqFileExist:
    print("Corresponding file containing all gap-filled sequences: " + gapfillSeqFile + "\n")


#----------------------------------------------------
# Summary Output
#----------------------------------------------------
try:
    print("\n------------------------------------------------------------------------------------------------------------------------\n")
    print("MTG-Link " + main.module)
    print("------------")

    # Go in the 'outDir' directory.
    try:
        os.chdir(main.outDir)
    except OSError as err:
        print("File 'mtglink.py': Something wrong with specified directory 'outDir'. \nOSError-{}".format(err))
        sys.exit(1)

    gfa_output = gfapy.Gfa.from_file(str(outputGFAFile))

    # Total initials gaps.
    totalGapsList = []
    for g_line in gfa.gaps:
        gap_start = str(g_line.sid1) +"_"+ str(g_line.sid2) 
        totalGapsList.append(gap_start)
    totalGapsNumber = len(totalGapsList)
    print("\nAttempt to gap-fill {} gaps \n".format(totalGapsNumber))

    # Gap(s) not gap-filled.
    gapsNotFilledList = []
    for g_line in gfa_output.gaps:
        gap_end = str(g_line.sid1) +"_"+ str(g_line.sid2) 
        gapsNotFilledList.append(gap_end)
        print("The gap {} was not successfully gap-filled".format(gap_end))

    gapsFilledNumber = len(totalGapsList) - len(gapsNotFilledList)
    print("\nIn total, {} gaps were successfully gap-filled:".format(str(gapsFilledNumber)))

    # Gaps gap-filled.
    if gapfillSeqFileExist:
        gapsNamesList = []
        if (gapfillSeqFile) is not None:
            try:
                with open(gapfillSeqFile, "r") as gapfilled:
                    for record in SeqIO.parse(gapfilled, "fasta"):
                        gapName_wo_rightSign = re.split('\+_|\-_', str(record.id))[0]
                        r2 = re.findall(r"\+_|\-_", str(record.id))[0]
                        gapName_rightSign = re.split("_", str(r2))[0]
                        gapName = gapName_wo_rightSign + gapName_rightSign

                        # For a new gap.
                        if gapName not in gapsNamesList:
                            gapsNamesList.append(gapName)
                            print("\t* " + gapName)

                        # For all gaps.
                        orientation = (re.split('\+_|\-_', str(record.id))[1]).split("_")[0]
                        length = str(record.description).split(' len.')[1].split('_qual.')[0]
                        quality = str(record.description).split('_qual.')[1]
                        print("\t\t* " + orientation + "\t" + length + " bp\t" + quality)
            
            except IOError as err:
                print("File 'mtglink.py', step 'Summary Output': Unable to open the file {}. \nIOError-{}".format(str(gapfillSeqFile), err))
                sys.exit(1)
            
    print("\n")

except Exception as e:
    print("File 'mtglink.py': Something wrong with the 'Summary Output' step.")
    print("Exception-{}".format(e))
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


#----------------------------------------------------
# Remove raw files from 'Local Assembly' step
#----------------------------------------------------
try:
    ## Local assembly with the DBG algorithm.
    if main.module == "DBG":

        # Go in the 'assemblyDir' directory.
        try:
            os.chdir(main.assemblyDir)
        except OSError as err:
            print("File 'mtglink.py': Something wrong with specified directory 'assemblyDir'. \nOSError-{}".format(err))
            sys.exit(1)
        
        # Remove the raw files.
        subprocess.run("rm -f *.h5", shell=True)
        subprocess.run("rm -f *.vcf", shell=True)
        subprocess.run("rm -f *.insertions.fasta", shell=True)

except Exception as e:
    print("File 'mtglink.py': Something wrong with removing the raw files from the 'Local Assembly' step.")
    print("Exception-{}".format(e))
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)

