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

"""Module 'mtglink.py': GFA processing and multi-threading

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
from main import gfaFile, gfa_name, line_gfa, outDir, subsamplingDir, assemblyDir, evalDir, module
from Pipeline import gapfilling
from helpers import update_gfa_with_solution


#----------------------------------------------------
# GFA pre-processing
#----------------------------------------------------
try:
    # Open the input GFA file.
    gfa = gfapy.Gfa.from_file(gfaFile)
    # Create the output GFA file.
    out_gfaFile = str(gfa_name).split('.gfa')[0] + "_mtglink.gfa"

    #----------------------------------------------------
    # GFA output: case no gap
    #----------------------------------------------------
    # If no gap, rewrite all the lines into GFA output.
    if len(gfa.gaps) == 0:
        with open(out_gfaFile, "w") as f:
            out_gfa = gfapy.Gfa()
            for line in gfa.lines:
                out_gfa.add_line(str(line))
            out_gfa.to_file(out_gfaFile)

    #----------------------------------------------------   
    # Fill the gaps
    #----------------------------------------------------
    # If gap, rewrite the H and S lines into GFA output.
    if line_gfa == "":
        with open(out_gfaFile, "w") as f:
            out_gfa = gfapy.Gfa()
            out_gfa.add_line("H\tVN:Z:2.0")
            for line in gfa.segments:
                out_gfa.add_line(str(line))
            out_gfa.to_file(out_gfaFile)
        
    gaps = []
    gaps_label = []
    # If '-line' argument provided, start analysis from this line in GFA file input.
    if line_gfa != "":
        for _gap_ in gfa.gaps[(line_gfa - (len(gfa.segments)+2)):]:
            _gap_ = str(_gap_)
            gaps.append(_gap_)
    else:
        # Convert Gfapy gap line to a string to be able to use it with multiprocessing.
        for _gap_ in gfa.gaps:
            _gap_ = str(_gap_)
            gaps.append(_gap_)

    p = Pool()

    with open("{}.union.sum".format(gfa_name), "w") as union_sum:
        legend = ["Gap_ID", "Left_scaffold", "Right_scaffold", "Gap_size", "Chunk_size", "Nb_barcodes", "Nb_reads"]
        union_sum.write('\t'.join(j for j in legend))

        for union_summary, output_for_gfa in p.map(gapfilling, gaps):
            # Write all union_summary (obtained for each gap) from 'gapfilling' into the 'union_sum' file.
            union_sum.write("\n" + '\t'.join(str(i) for i in union_summary))

            # Output the 'output_for_gfa' results (obtained for each gap) from 'gapfilling' in the output GFA file.
            print("\nCreating the output GFA file...")
            ## Solution found for the current gap
            if len(output_for_gfa[0]) > 1:          
                for output in output_for_gfa:
                    gapfillFile = update_gfa_with_solution(outDir, gfa_name, output, out_gfaFile)
                    success = True
            ## No solution found for the current gap
            else:                                   
                out_gfa = gfapy.Gfa.from_file(out_gfaFile)
                out_gfa.add_line(output_for_gfa[0][0])
                out_gfa.to_file(out_gfaFile)
                success = False

        p.close()

    # # Remove the raw files obtained from MindTheGap.
    # os.chdir(assemblyDir)
    # subprocess.run("rm -f *.h5", shell=True)
    # subprocess.run("rm -f *.vcf", shell=True)

except Exception as e:
    print("\nException-")
    print(e)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)


print("\nThe results from MTG-Link are saved in: " + outDir)
print("The results from the Read Subsampling step are saved in " + subsamplingDir)
print("The results from the Local Assembly step are saved in " + assemblyDir)
print("The results from the Qualitative Evaluation step are saved in " + evalDir)

print("\nSummary of the union: " + gfa_name + ".union.sum")
print("GFA output file: " + out_gfaFile)
if success == True:
    print("Corresponding file containing all gap-filled sequences: " + gapfillFile + "\n")


#----------------------------------------------------
#Summary output
#----------------------------------------------------
try:
    print("\n------------------------------------------------------------------------------------------------------------------------\n")
    print("MTG-Link " + module)
    print("------------")

    gfa_output = gfapy.Gfa.from_file(outDir +"/"+ str(out_gfaFile))

    # Total initials gaps.
    total_gaps = []
    for g_line in gfa.gaps:
        gap_start = str(g_line.sid1) +"_"+ str(g_line.sid2) 
        total_gaps.append(gap_start)
    nb_total_gaps = len(total_gaps)
    print("\nAttempt to gap-fill {} gaps \n".format(nb_total_gaps))

    # Gap(s) not gap-filled.
    no_gapfill = []
    for g_line in gfa_output.gaps:
        gap_end = str(g_line.sid1) +"_"+ str(g_line.sid2) 
        no_gapfill.append(gap_end)
        print("The gap {} was not successfully gap-filled".format(gap_end))

    nb_gapfill = len(total_gaps) - len(no_gapfill)
    print("In total, {} gaps were successfully gap-filled:".format(str(nb_gapfill)))

    # Gaps gap-filled.
    out_fasta_file = outDir +"/"+ gapfillFile
    gap_names = []
    if (out_fasta_file) is not None:
        with open(out_fasta_file, "r") as gapfilled:
            for record in SeqIO.parse(gapfilled, "fasta"):
                gap_name = str(record.id).split('_')[0]

                # For a new gap.
                if gap_name not in gap_names:
                    gap_names.append(gap_name)
                    print("\t* " + gap_name)

                # For all gaps.
                orientation = str(record.id).split('_')[-1]
                length = str(record.description).split('_ len_')[1].split('_qual_')[0]
                quality = str(record.description).split('_qual_')[1]
                print("\t\t* " + orientation + "\t" + length + " bp\t" + quality)
            
    print("\n")

except Exception as e:
    print("\nException-")
    print(e)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)

