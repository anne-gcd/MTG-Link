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

from __future__ import print_function
import os
import sys
import re
import csv
import argparse
import subprocess
from Bio import Align


# PairwiseAligner object.
aligner = Align.PairwiseAligner()
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.query_end_gap_score = 0
aligner.internal_open_gap_score = -1
aligner.internal_extend_gap_score = -0.5
aligner.target_end_open_gap_score = -1
aligner.target_end_extend_gap_score = -0.5


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="stats_alignment.py", usage="%(prog)s -qry <query_sequences_file> -ref <reference_sequence> -ext <extension_size> -p <output_file_prefix> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Statistics about the Nucmer alignments between a Query and a Reference
                                Note: there are k-mer flanking regions on the edges of the assembled sequence (Query) 
                                Note: the assembled sequences are those for the extended target (so they include the '-ext' bp flanking regions)
                                '''))

#TODO: Modif lines 55 and 63-64
parser.add_argument("-qry", "--query", action="store", help="File containing the query sequences (format: 'xxx.insertions_filtered.fasta')", required=True)
parser.add_argument("-ref", "--reference", action="store", help="File containing either the reference sequence or the target flanking contigs sequences (format: 'xxx.fasta')", required=True)
parser.add_argument("-ext", "--ext", action="store", type=int, help="Size of the extension of the target on both sides (bp); determine start/end of local assembly", required=True)
parser.add_argument("-p", "--prefix", action="store", help="Prefix of the output files", required=True)
parser.add_argument("-out", "--outDir", action="store", help="Output directory for saving results", required=True)

args = parser.parse_args()

if re.match('^.*.insertions_filtered.fasta$', args.query) is None:
    parser.error("\nWarning: Qualitative evaluation _ The suffix of the query sequences file should be: '.insertions_filtered.fasta'")

if re.match('^.*.fasta$', args.reference) is None:
    parser.error("\nWarning: Qualitative evaluation _ The suffix of the reference sequence file should be: '.fasta'")


#----------------------------------------------------
# Input files
#----------------------------------------------------
# Query file: assembled sequences file.
qry_file = os.path.abspath(args.query)
if not os.path.exists(args.query):
    parser.error("\nWarning: Qualitative evaluation _ The path of the query file doesn't exist")

# Reference file: containing either the reference sequence or the flanking contigs sequences.
ref_file = os.path.abspath(args.reference)
if not os.path.exists(ref_file):
    parser.error("\nWarning: Qualitative evaluation _ The path of the reference file doesn't exist")


#----------------------------------------------------
# Directory for saving results
#----------------------------------------------------
cwd = os.getcwd()

# Create the directory 'outDir', where the results will be saved.
if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)
try:
    os.chdir(args.outDir)
except OSError:
    print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()


# #----------------------------------------------------
# # Ref = reference sequence of simulated gap/target
# #----------------------------------------------------
# if not re.match('^.*.contigs.fasta$', args.reference):

#     #----------------------------------------------------
#     # NUCmer alignments
#     #----------------------------------------------------
#     try: 
#         # Run NUCmer to obtain alignments of the reference sequence against the query's sequences.
#         prefix = args.prefix + ".ref_qry"

#         log_file = str(prefix) + ".log"
#         with open(log_file, "a") as log:
#             log.write("Query file: " + str(qry_file) + "\n")
#             log.write("Reference file: " + str(ref_file) + "\n")
#             log.write("The results are saved in " + outDir)

#         nucmerLog = "{}_nucmer_ref_qry.log".format(args.prefix)
#         delta_file = prefix + ".delta"
#         coords_file = prefix + ".coords.unsorted" 

#         # Keep only alignments with >90% Id. ('-I90'). 
#         #nucmer_command = ["nucmer", "--maxmatch", "-p", prefix, ref_file, qry_file]
#         nucmer_command = ["nucmer", "-p", prefix, ref_file, qry_file]
#         coords_command = ["show-coords", "-rcdlT", "-I90", delta_file]

#         with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
#             subprocess.run(nucmer_command, stderr=log)
#             subprocess.run(coords_command, stdout=coords, stderr=log)

#         # Sort the 'xxx.coords.unsorted' file for further analysis.
#         coords_sorted_file = prefix + ".coords"
#         sort_command = ["sort", "-n", coords_file]
#         with open(coords_sorted_file, "w") as coords_sorted:
#             subprocess.run(sort_command, stdout=coords_sorted)

#     except Exception as e:
#         print("\nFile 'stats_alignment.py': Something wrong with the NUCmer alignments, when ref = reference sequence")
#         print("Exception-")
#         print(e)
#         sys.exit(1)

#     #----------------------------------------------------
#     # Parameters of the local assembly step
#     #----------------------------------------------------
#     try:
#         # Local assembly performed with the DBG algorithm.
#         if ".k" in qry_file.split('/')[-1]:
#             gap_size = qry_file.split('.bxu')[0].split('.')[-5]
#             g = int("".join(list(gap_size)[1:]))
#             if g == 0:
#                 g = "NA"
#             flank_size = qry_file.split('.bxu')[0].split('.')[-4]
#             flank = int(flank_size.split('flank')[1])
#             barcodes_occ = qry_file.split('.bxu')[0].split('.')[-3]
#             occ = int(barcodes_occ.split('occ')[1])
#             kmer_size = qry_file.split('.bxu')[0].split('.')[-2]
#             k = int("".join(list(kmer_size)[1:]))
#             abundance_min = qry_file.split('.bxu')[0].split('.')[-1]
#             a = int("".join(list(abundance_min)[1:]))
#             qry_id = qry_file.split('/')[-1].split('.')[2]

#         # Local assembly performed with the IRO algorithm.
#         if ".dmax" in qry_file.split('/')[-1]:
#             gap_size = qry_file.split('.bxu')[0].split('.')[-7]
#             g = int("".join(list(gap_size)[1:]))
#             if g == 0:
#                 g = "NA"
#             flank_size = qry_file.split('.bxu')[0].split('.')[-6]
#             flank = int(flank_size.split('flank')[1])
#             barcodes_occ = qry_file.split('.bxu')[0].split('.')[-5]
#             occ = int(barcodes_occ.split('occ')[1])
#             seed_size = qry_file.split('.bxu')[0].split('.')[-4]
#             s = int("".join(list(seed_size)[1:]))
#             min_overlap = qry_file.split('.bxu')[0].split('.')[-3]
#             o = int("".join(list(min_overlap)[1:]))
#             abundance_min = qry_file.split('.bxu')[0].split('.')[-2]
#             a = str("".join(list(abundance_min)[1:]))
#             dmax = qry_file.split('.bxu')[0].split('.')[-1]
#             d = int("".join(list(dmax)[4:]))
#             qry_id = qry_file.split('/')[-1].split('.')[2]

#     except Exception as e:
#         print("\nFile 'stats_alignment.py': Something wrong with getting the parameters of the local assembly step, when ref = reference sequence")
#         print("Exception-")
#         print(e)
#         sys.exit(1)

#     #----------------------------------------------------
#     # Quality estimation
#     #----------------------------------------------------
#     try:
#         # Output stats file of NUCmer alignment (Query vs Ref).
#         ref_qry_output = outDir + "/" + args.prefix + ".nucmerAlignments.stats.unsorted"

#         ## Local assembly performed with the DBG algorithm
#         if ".k" in qry_file.split('/')[-1]:
#             stats_legend = ["Target", "Len_target", "Flank", "Barc_occ", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
#                             "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]

#         ## Local assembly performed with the IRO algorithm
#         if ".dmax" in qry_file.split('/')[-1]:
#             stats_legend = ["Target", "Len_target", "Flank", "Barc_occ", "s", "o", "a", "dmax", "Len_Q", "Ref", "Len_R", \
#                             "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]
        
#         # Get output values from NUCmer.
#         reader = csv.DictReader(open(coords_sorted_file), \
#                                     fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
#                                     delimiter='\t')

#         rows = list(reader)
#         for row in rows[3:]:
#                 len_q = row["LEN_Q"]
#                 ref = row["TAG_1"]
#                 len_r = row["LEN_R"]
#                 start_r = row["S1"]
#                 end_r = row["E1"]
#                 start_q = row["S2"]
#                 end_q = row["E2"]
#                 len_align_r = row["LEN_1"]
#                 len_align_q = row["LEN_2"]
#                 identity = row["%_IDY"]
#                 cov_r = row["COV_R"]
#                 cov_q = row["COV_Q"]
#                 frame_r = row["FRM_R"]
#                 frame_q = row["FRM_Q"]

#                 # Local assembly performed with the DBG algorithm.
#                 if ".k" in qry_file.split('/')[-1]:
#                     solution = str(row["TAG_2"]).split('_sol_')[1]
#                     if "bkpt1" in str(row["TAG_2"]):
#                         strand = "fwd"
#                     else:
#                         strand = "rev"

#                 # Estimate quality of assembled sequence (Query).
#                 ref_len = int(len_r)
#                 qry_len = int(len_q) - 2*args.ext
#                 error_10_perc = int(0.1 * ref_len)
#                 error_50_perc = int(0.5 * ref_len)

#                 # Length of query sequence is equal +-10% of ref length.
#                 if qry_len in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
#                     ## The assembled seq matches to the whole ref seq
#                     if int(len_align_q) == ref_len:
#                         quality_rq = 'A'
#                     ## The assembled seq matches to the ref seq +-10% of ref length
#                     elif int(len_align_q) in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
#                         quality_rq = 'B'
#                     ## The assembled seq matches to the ref seq, but not along all their length (>= 50% of their length align)
#                     elif int(len_align_q) in range((ref_len - error_50_perc), (ref_len + error_50_perc)):
#                         quality_rq = 'C'
#                     else:
#                         quality_rq = 'D'

#                 else:
#                     quality_rq = 'D'

#                 # Write stats results in output file.
#                 ## Local assembly performed with the DBG algorithm
#                 if ".k" in qry_file.split('/')[-1]:
#                     stats = [qry_id, g, flank, occ, k, a, strand, solution, len_q, ref, len_r, \
#                         start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

#                 ## Local assembly performed with the IRO algorithm
#                 if ".dmax" in qry_file.split('/')[-1]:
#                     stats = [qry_id, g, flank, occ, s, o, a, d, len_q, ref, len_r, \
#                         start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

#                 ## Write in output file
#                 if os.path.exists(ref_qry_output):
#                     with open(ref_qry_output, "a") as output:
#                         output.write('\n' + '\t'.join(str(i) for i in stats))
#                 else:
#                     with open(ref_qry_output, "a") as output:
#                         output.write('\t'.join(j for j in stats_legend))
#                         output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))

#         # Sort the 'xxx.nucmerAlignments.stats.unsorted' file for further analysis.
#         ref_qry_sorted = outDir + "/" + args.prefix + ".nucmerAlignments.stats"
#         order_command = ["sort", "-k7,8", "-k10", "-k12,13n", "-r", ref_qry_output]
#         with open(ref_qry_sorted, "w") as r_sorted:
#             subprocess.run(order_command, stdout=r_sorted)

#     except Exception as e:
#         print("\nFile 'stats_alignment.py': Something wrong with the quality estimation, when ref = reference sequence")
#         print("Exception-")
#         print(e)
#         sys.exit(1)
        
#     #----------------------------------------------------
#     # Alignment in multiple chunks/flanks (for DBG)
#     #----------------------------------------------------
#     try:      
#         # If alignment in multiple chunks/flanks, calculate the appropriate quality score.
#         ## Local assembly performed with the DBG algorithm
#         if ".k" in qry_file.split('/')[-1]:
#             with open(ref_qry_sorted, "r") as r:
#                 r.seek(0)
#                 reader = csv.DictReader(r, fieldnames=("Target", "Len_target", "Flank", "Barc_occ", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
#                                         "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
#                                         delimiter='\t')

#                 rows = list(reader)
#                 for i in range(1, len(rows)):
#                     if ("fwd" in rows[i]["Strand"]) or ("rev" in rows[i]["Strand"]):
#                         qry_id = rows[i]["Target"]
#                         g = rows[i]["Len_target"]
#                         c = rows[i]["Flank"]
#                         f = rows[i]["Barc_occ"]
#                         k = rows[i]["k"]
#                         a = rows[i]["a"]
#                         strand = rows[i]["Strand"]
#                         solution = rows[i]["Solution"]
#                         len_q = rows[i]["Len_Q"]
#                         ref = rows[i]["Ref"]
#                         len_r = rows[i]["Len_R"]
#                         frame_r = rows[i]["Frame_R"]
#                         frame_q = rows[i]["Frame_Q"]
#                         ref_len = int(len_r)
#                         qry_len = int(len_q) - 2*args.ext

#                         if (strand != rows[i-1]["Strand"]) or (solution != rows[i-1]["Solution"]):
#                             lack_ref = 0
#                             lack_qry = 0
#                             start_r = int(rows[i]["Start_ref"])
#                             lack_ref += start_r - 1
#                             if frame_q == "1":
#                                 start_q = int(rows[i]["Start_qry"])
#                                 lack_qry += start_q - (args.ext+1)
#                             elif frame_q == "-1":
#                                 end_q = int(rows[i]["Start_qry"])
#                                 lack_qry += qry_len - (end_q - args.ext)

#                         if (strand == rows[i-1]["Strand"]) and (solution == rows[i-1]["Solution"]):
#                             if int(rows[i]["Start_ref"]) > int(rows[i-1]["End_ref"]):
#                                 lack_ref += int(rows[i]["Start_ref"]) - int(rows[i-1]["End_ref"])
#                             if frame_q == "1":
#                                 if int(rows[i]["Start_qry"]) > int(rows[i-1]["End_qry"]):
#                                     lack_qry += int(rows[i]["Start_qry"]) - int(rows[i-1]["End_qry"])
#                             elif frame_q == "-1":
#                                 if int(rows[i]["Start_qry"]) < int(rows[i-1]["End_qry"]):
#                                     lack_qry += int(rows[i-1]["End_qry"]) - int(rows[i]["Start_qry"])

#                         if (i == len(rows)-1) or ((strand != rows[i+1]["Strand"]) or (solution != rows[i+1]["Solution"])):
#                             end_r = int(rows[i]["End_ref"])
#                             lack_ref += ref_len - end_r
#                             if frame_q == "1":
#                                 end_q = int(rows[i]["End_qry"])
#                                 lack_qry += qry_len - (end_q - args.ext)
#                             elif frame_q == "-1":
#                                 start_q = int(rows[i]["End_qry"])
#                                 lack_qry += start_q - (args.ext+1)

#                             len_align_r = ref_len - lack_ref
#                             len_align_q = qry_len - lack_qry

#                             identity = "/"
#                             cov_r = "/"
#                             cov_q = "/"

#                             # Assign a quality score.
#                             ## The assembled seq matches to the whole ref seq
#                             if int(len_align_q) == ref_len:
#                                 quality_rq = 'A'
#                             ## The assembled seq matches to the ref seq +-10% of ref length
#                             elif int(len_align_q) in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
#                                 quality_rq = 'B'
#                             ## The assembled seq matches to the ref seq, but not along all their length (>= 50% of their length align)
#                             elif int(len_align_q) in range((ref_len - error_50_perc), (ref_len + error_50_perc)):
#                                 quality_rq = 'C'
#                             else:
#                                 quality_rq = 'D'

#                             # Write stats results in output file.
#                             ## Local assembly performed with the DBG algorithm
#                             if ".k" in qry_file.split('/')[-1]:
#                                 stats = [qry_id, g, flank, occ, k, a, strand, solution, len_q, ref, len_r, \
#                                     start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

#                             ## Local assembly performed with the IRO algorithm
#                             if ".dmax" in qry_file.split('/')[-1]:
#                                 stats = [qry_id, g, flank, occ, s, o, a, d, len_q, ref, len_r, \
#                                     start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

#                             with open(ref_qry_sorted, "a") as output_ref:
#                                 output_ref.write('\n' + '\t'.join(str(i) for i in stats))

#     except Exception as e:
#         print("\nFile 'stats_alignment.py': Something wrong with the quality estimation for alignment in multiple chunks/flanks (DBG), when ref = reference sequence")
#         print("Exception-")
#         print(e)
#         sys.exit(1)
        

#----------------------------------------------------
# Ref = gap/target flanking contigs' sequences
#----------------------------------------------------
# elif re.match('^.*.contigs.fasta$', args.reference):
if re.match('^.*.contigs.fasta$', args.reference):

    #----------------------------------------------------
    # NUCmer alignments
    #----------------------------------------------------
    try:
        # Run NUCmer to obtain alignments of extension portions (-ext) (of the flanking contigs) against the query's sequences.
        prefix = args.prefix + ".ref_qry"

        log_file = str(prefix) + ".log"
        with open(log_file, "a") as log:
            log.write("Query file: " + str(qry_file) + "\n")
            log.write("Reference file: " + str(ref_file) + "\n")
            log.write("The results are saved in " + outDir)

        nucmerLog = "{}_nucmer_ref_qry.log".format(args.prefix)
        delta_file = prefix + ".delta"
        coords_file = prefix + ".coords.unsorted"

        # Keep only alignments with >90% Id. ('-I90'). 
        nucmer_command = ["nucmer", "-p", prefix, ref_file, qry_file]
        coords_command = ["show-coords", "-rcdlT", "-I90", delta_file]

        with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
            subprocess.run(nucmer_command, stderr=log)
            subprocess.run(coords_command, stdout=coords, stderr=log)

        # Sort the 'xxx.coords.unsorted' file for further analysis.
        coords_sorted_file = prefix + ".coords"
        sort_command = ["sort", "-n", coords_file]
        with open(coords_sorted_file, "w") as coords_sorted:
            subprocess.run(sort_command, stdout=coords_sorted)

    except Exception as e:
        print("\nFile 'stats_alignment.py': Something wrong with the NUCmer alignments, when ref = flanking contigs' sequences")
        print("Exception-")
        print(e)
        sys.exit(1)

    #----------------------------------------------------
    # Parameters of the local assembly step
    #----------------------------------------------------
    try:
        # Local assembly performed with the DBG algorithm.
        if ".k" in qry_file.split('/')[-1]:
            # file name of the form: test.gfa.27358_0-171884-L+_27358_172884-344768-R+.g1000.flank10000.occ2.k61.a3.bxu.insertions_filtered.fasta
            gap_size = qry_file.split('.bxu')[0].split('.')[-5]
            g = int("".join(list(gap_size)[1:]))
            if g == 0:
                g = "NA"
            flank_size = qry_file.split('.bxu')[0].split('.')[-4]
            flank = int(flank_size.split('flank')[1])
            barcodes_occ = qry_file.split('.bxu')[0].split('.')[-3]
            occ = int(barcodes_occ.split('occ')[1])
            kmer_size = qry_file.split('.bxu')[0].split('.')[-2]
            k = int("".join(list(kmer_size)[1:]))
            abundance_min = qry_file.split('.bxu')[0].split('.')[-1]
            a = int("".join(list(abundance_min)[1:]))
            qry_id = qry_file.split('/')[-1].split('.gfa.')[-1].split('.')[0]

        # Local assembly performed with the IRO algorithm.
        if ".dmax" in qry_file.split('/')[-1]:
            gap_size = qry_file.split('.bxu')[0].split('.')[-7]
            g = int("".join(list(gap_size)[1:]))
            if g == 0:
                g = "NA"
            flank_size = qry_file.split('.bxu')[0].split('.')[-6]
            flank = int(flank_size.split('flank')[1])
            barcodes_occ = qry_file.split('.bxu')[0].split('.')[-5]
            occ = int(barcodes_occ.split('occ')[1])
            seed_size = qry_file.split('.bxu')[0].split('.')[-4]
            s = int("".join(list(seed_size)[1:]))
            min_overlap = qry_file.split('.bxu')[0].split('.')[-3]
            o = int("".join(list(min_overlap)[1:]))
            abundance_min = qry_file.split('.bxu')[0].split('.')[-2]
            a = str("".join(list(abundance_min)[1:]))
            dmax = qry_file.split('.bxu')[0].split('.')[-1]
            d = int("".join(list(dmax)[4:]))
            qry_id = qry_file.split('/')[-1].split('.gfa.')[-1].split('.')[0]

    except Exception as e:
        print("\nFile 'stats_alignment.py': Something wrong with getting the parameters of the local assembly step, when ref = flanking contigs' sequences")
        print("Exception-")
        print(e)
        sys.exit(1)

    #----------------------------------------------------
    # Quality estimation
    #----------------------------------------------------
    try:
        # Output stats file of NUCmer alignment (Query vs Ref).
        ref_qry_output = outDir + "/" + args.prefix + ".nucmerAlignments.stats.unsorted"

        ## Local assembly performed with the DBG algorithm
        if ".k" in qry_file.split('/')[-1]:
            stats_legend = ["Target", "Len_target", "Chunk", "Barc_occ", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                            "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]

        ## Local assembly performed with the IRO algorithm
        if ".dmax" in qry_file.split('/')[-1]:
            stats_legend = ["Target", "Len_target", "Chunk", "Barc_occ", "s", "o", "a", "dmax", "Len_Q", "Ref", "Len_R", \
                            "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]
        
        ## Write legend in output file.
        with open(ref_qry_output, "a") as output:
            output.write('\t'.join(j for j in stats_legend))
        
        # Get output values from NUCmer.
        reader = csv.DictReader(open(coords_sorted_file), \
                                    fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                    delimiter='\t')

        rows = list(reader)
        for row in rows[3:]:
            if row["TAG_1"].split("_region")[0] in str(qry_id):
                len_q = row["LEN_Q"]
                ref = row["TAG_1"].split("_region")[0]
                len_r = row["LEN_R"]
                start_r = row["S1"]
                end_r = row["E1"]
                start_q = row["S2"]
                end_q = row["E2"]
                len_align_r = row["LEN_1"]
                len_align_q = row["LEN_2"]
                identity = row["%_IDY"]
                cov_r = row["COV_R"]
                cov_q = row["COV_Q"]
                frame_r = row["FRM_R"]
                frame_q = row["FRM_Q"]

                # Local assembly performed with the DBG algorithm.
                if ".k" in qry_file.split('/')[-1]:
                    solution = str(row["TAG_2"]).split('_sol_')[1]
                    if "bkpt1" in str(row["TAG_2"]):
                        strand = "fwd"
                    else:
                        strand = "rev"

                # Estimate quality of assembled sequence (Query).
                left_wo_sign = re.split('\+_|\-_', str(qry_id))[0]
                r1 = re.findall(r"\+_|\-_",str(qry_id))[0]
                left_sign = re.split("_", str(r1))[0]
                left = left_wo_sign + left_sign
                left_scaffold = left[:-1]
                right = re.split('\+_|\-_', str(qry_id))[1]
                right_scaffold = right[:-1]
                error_10_perc = int(0.1 * args.ext)
                error_50_perc = int(0.5 * args.ext)

                # Local assembly performed with the DBG algorithm.
                if ".k" in qry_file.split('/')[-1]:
    
                    # Ref = Left scaffold.
                    if ref == left_scaffold:                    
                        ## Extension of qry match perfectly as expected to ref
                        if ('+' in left and ((strand == "fwd" and int(start_q) == 1 and int(end_q) == args.ext) or (strand == "rev" and int(start_q) == int(len_q) and int(end_q) == (int(len_q) - args.ext + 1)))) \
                            or ('-' in left and ((strand == "fwd" and int(start_q) == args.ext and int(end_q) == 1) or (strand == "rev" and int(start_q) == (int(len_q) - args.ext + 1) and int(end_q) == int(len_q)))):
                            quality_rq = 'A'
                        ## Extension of qry almost match as expected to ref (+-10% of extension size) 
                        elif ('+' in left and ((strand == "fwd" and int(start_q) in range(1, (1 + error_10_perc+1)) and int(end_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1))) or (strand == "rev" and int(start_q) in range((int(len_q) - error_10_perc), (int(len_q)+1)) and int(end_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1))))) \
                            or ('-' in left and ((strand == "fwd" and int(start_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1)) and int(end_q) in range(1, (1 + error_10_perc+1))) or (strand == "rev" and int(start_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1)) and int(end_q) in range((int(len_q) - error_10_perc), (int(len_q)+1))))):
                            quality_rq = 'B'
                        ## Extension of qry almost match as expected to ref (+-50% of extension size) 
                        elif ('+' in left and ((strand == "fwd" and int(start_q) in range(1, (1 + error_50_perc+1)) and int(end_q) in range((args.ext - error_50_perc), (args.ext + error_50_perc+1))) or (strand == "rev" and int(start_q) in range((int(len_q) - error_50_perc), (int(len_q)+1)) and int(end_q) in range((int(len_q)-args.ext+1 - error_50_perc), (int(len_q)-args.ext+1 + error_50_perc+1))))) \
                            or ('-' in left and ((strand == "fwd" and int(start_q) in range((args.ext - error_50_perc), (args.ext + error_50_perc+1)) and int(end_q) in range(1, (1 + error_50_perc+1))) or (strand == "rev" and int(start_q) in range((int(len_q)-args.ext+1 - error_50_perc), (int(len_q)-args.ext+1 + error_50_perc+1)) and int(end_q) in range((int(len_q) - error_50_perc), (int(len_q)+1))))):
                            quality_rq = 'C'
                        else:
                            quality_rq = 'D'

                    # Ref = Right scaffold.
                    elif ref == right_scaffold:
                        ## Extension of qry match perfectly as expected to ref
                        if ('+' in right and ((strand == "fwd" and int(start_q) == (int(len_q) - args.ext + 1) and int(end_q) == int(len_q)) or (strand == "rev" and int(start_q) == args.ext and int(end_q) == 1))) \
                            or ('-' in right and ((strand == "fwd" and int(start_q) == int(len_q) and int(end_q) == (int(len_q) - args.ext +1)) or (strand == "rev" and int(start_q) == 1 and int(end_q) == args.ext))):
                            quality_rq = 'A'
                        ## Extension of qry almost match as expected to ref (+-10% of extension size) 
                        elif ('+' in right and ((strand == "fwd" and int(start_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1)) and int(end_q) in range((int(len_q) - error_10_perc), (int(len_q)+1))) or (strand == "rev" and int(start_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1)) and int(end_q) in range(1, (1 + error_10_perc+1))))) \
                            or ('-' in right and ((strand == "fwd" and int(start_q) in range((int(len_q) - error_10_perc), (int(len_q)+1)) and int(end_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1))) or (strand == "rev" and int(start_q) in range(1, (1 + error_10_perc+1)) and int(end_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1))))):
                            quality_rq = 'B'
                        ## Extension of qry almost match as expected to ref (+-50% of extension size) 
                        elif ('+' in right and ((strand == "fwd" and int(start_q) in range((int(len_q)-args.ext+1 - error_50_perc), (int(len_q)-args.ext+1 + error_50_perc+1)) and int(end_q) in range((int(len_q) - error_50_perc), (int(len_q)+1))) or (strand == "rev" and int(start_q) in range((args.ext - error_50_perc), (args.ext + error_50_perc+1)) and int(end_q) in range(1, (1 + error_50_perc+1))))) \
                            or ('-' in right and ((strand == "fwd" and int(start_q) in range((int(len_q) - error_50_perc), (int(len_q)+1)) and int(end_q) in range((int(len_q)-args.ext+1 - error_50_perc), (int(len_q)-args.ext+1 + error_50_perc+1))) or (strand == "rev" and int(start_q) in range(1, (1 + error_50_perc+1)) and int(end_q) in range((args.ext - error_50_perc), (args.ext + error_50_perc+1))))):
                            quality_rq = 'C'
                        else:
                            quality_rq = 'D'

                # Local assembly performed with the IRO algorithm.
                if ".dmax" in qry_file.split('/')[-1]:
                    
                    # Ref = Left scaffold.
                    if ref == left_scaffold:
                        ## Extension of qry match perfectly as expected to ref
                        if ('+' in left and int(start_q) == 1 and int(end_q) == args.ext) \
                            or ('-' in left and int(start_q) == args.ext and int(end_q) == 1):
                            quality_rq = 'A'
                        ## Extension of qry almost match as expected to ref (+-10% of extension size) 
                        elif ('+' in left and int(start_q) in range(1, (1 + error_10_perc + 1)) and int(end_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc + 1))) \
                            or ('-' in left and int(start_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc + 1)) and int(end_q) in range(1, (1 + error_10_perc + 1))):
                            quality_rq = 'B'
                        ## Extension of qry almost match as expected to ref (+-50% of extension size) 
                        elif ('+' in left and int(start_q) in range(1, (1 + error_50_perc + 1)) and int(end_q) in range((args.ext - error_50_perc), (args.ext + error_50_perc + 1))) \
                            or ('-' in left and int(start_q) in range((args.ext - error_50_perc), (args.ext + error_50_perc + 1)) and int(end_q) in range(1, (1 + error_50_perc + 1))):
                            quality_rq = 'C'
                        else:
                            quality_rq = 'D'

                    # Ref = Right scaffold.
                    elif ref == right_scaffold:
                        ## Extension of qry match perfectly as expected to ref
                        if ('+' in right and int(start_q) == (int(len_q) - args.ext + 1) and int(end_q) == (int(len_q))) \
                            or ('-' in right and int(start_q) == (int(len_q)) and int(end_q) == (int(len_q) - args.ext + 1)):
                            quality_rq = 'A'
                        ## Extension of qry almost match as expected to ref (+-10% of extension size) 
                        elif ('+' in right and int(start_q) in range((int(len_q) - args.ext + 1 - error_10_perc), (int(len_q) - args.ext + 1 + error_10_perc + 1)) and int(end_q) in range((int(len_q) - error_10_perc), (int(len_q) + 1))) \
                            or ('-' in right and int(start_q) in range((int(len_q) - error_10_perc), (int(len_q) + 1)) and int(end_q) in range((int(len_q) - args.ext + 1 - error_10_perc), (int(len_q) - args.ext + 1 + error_10_perc + 1))):
                            quality_rq = 'B'
                        ## Extension of qry almost match as expected to ref (+-50% of extension size) 
                        elif ('+' in right and int(start_q) in range((int(len_q) - args.ext + 1 - error_50_perc), (int(len_q) - args.ext + 1 + error_50_perc + 1)) and int(end_q) in range((int(len_q) - error_50_perc), (int(len_q) + 1))) \
                            or ('-' in right and int(start_q) in range((int(len_q) - error_50_perc), (int(len_q) + 1)) and int(end_q) in range((int(len_q) - args.ext + 1 - error_50_perc), (int(len_q) - args.ext + 1 + error_50_perc + 1))):
                            quality_rq = 'C'
                        else:
                            quality_rq = 'D'

                # Write stats results in output file.
                ## Local assembly performed with the DBG algorithm
                if ".k" in qry_file.split('/')[-1]:
                    stats = [qry_id, g, flank, occ, k, a, strand, solution, len_q, ref, len_r, \
                        start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                ## Local assembly performed with the IRO algorithm
                if ".dmax" in qry_file.split('/')[-1]:
                    stats = [qry_id, g, flank, occ, s, o, a, d, len_q, ref, len_r, \
                        start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                ## Write in output file
                if os.path.exists(ref_qry_output):
                    with open(ref_qry_output, "a") as output:
                        output.write('\n' + '\t'.join(str(i) for i in stats))
                else:
                    with open(ref_qry_output, "a") as output:
                        output.write('\t'.join(j for j in stats_legend))
                        output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))

        # Sort the 'xxx.nucmerAlignments.stats.unsorted' file for further analysis.
        ref_qry_sorted = outDir + "/" + args.prefix + ".nucmerAlignments.stats"
        order_command = ["sort", "-k7,8", "-k10", "-k12,13n", "-r", ref_qry_output]
        with open(ref_qry_sorted, "w") as r_sorted:
            subprocess.run(order_command, stdout=r_sorted)

    except Exception as e:
        print("\nFile 'stats_alignment.py': Something wrong with the quality estimation, when ref = flanking contigs' sequences")
        print("Exception-")
        print(e)
        sys.exit(1)

    #----------------------------------------------------
    # Alignment in multiple chunks/flanks (for DBG)
    #----------------------------------------------------
    try:      
        # If alignment in multiple chunks/flanks, calculate the appropriate quality score.
        ## Local assembly performed with the DBG algorithm
        if ".k" in qry_file.split('/')[-1]:
            with open(ref_qry_sorted, "r") as r:
                r.seek(0)
                reader = csv.DictReader(r, fieldnames=("Target", "Len_target", "Flank", "Barc_occ", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
                                        delimiter='\t')

                rows = list(reader)
                for i in range(1, len(rows)):
                    if ("fwd" in rows[i]["Strand"]) or ("rev" in rows[i]["Strand"]):
                        qry_id = rows[i]["Target"]
                        g = rows[i]["Len_target"]
                        c = rows[i]["Flank"]
                        f = rows[i]["Barc_occ"]
                        k = rows[i]["k"]
                        a = rows[i]["a"]
                        strand = rows[i]["Strand"]
                        solution = rows[i]["Solution"]
                        len_q = rows[i]["Len_Q"]
                        ref = rows[i]["Ref"]
                        len_r = rows[i]["Len_R"]
                        frame_r = rows[i]["Frame_R"]
                        frame_q = rows[i]["Frame_Q"]
                        ref_len = int(len_r)
                        qry_len = int(len_q) - 2*args.ext

                        if (strand != rows[i-1]["Strand"] or solution != rows[i-1]["Solution"] or ref != rows[i-1]["Ref"]):
                            lack_ref = 0
                            lack_qry = 0
                            end_r = int(rows[i]["End_ref"])
                            lack_ref += ref_len - end_r
                            if frame_q == "1":
                                end_q = int(rows[i]["End_qry"])
                                lack_qry += qry_len - (end_q - args.ext)
                            elif frame_q == "-1":
                                start_q = int(rows[i]["End_qry"])
                                lack_qry += start_q - (args.ext+1)

                        if (strand == rows[i-1]["Strand"] and solution == rows[i-1]["Solution"] and ref == rows[i-1]["Ref"]):
                            if int(rows[i]["End_ref"]) < int(rows[i-1]["Start_ref"]):
                                lack_ref += (int(rows[i-1]["Start_ref"]) - int(rows[i]["End_ref"]))
                            if frame_q == "1":
                                if int(rows[i]["End_qry"]) < int(rows[i-1]["Start_qry"]):
                                    lack_qry += (int(rows[i-1]["Start_qry"]) - int(rows[i]["End_qry"]))
                            elif frame_q == "-1":
                                if int(rows[i]["End_qry"]) > int(rows[i-1]["Start_qry"]):
                                    lack_qry += (int(rows[i]["End_qry"]) - int(rows[i-1]["Start_qry"]))

                        if (i == len(rows)-1) or (strand != rows[i+1]["Strand"] or solution != rows[i+1]["Solution"] or ref != rows[i+1]["Ref"]):
                            start_r = int(rows[i]["Start_ref"])
                            lack_ref += start_r - 1
                            if frame_q == "1":
                                start_q = int(rows[i]["Start_qry"])
                                lack_qry += start_q - (args.ext+1)
                            elif frame_q == "-1":
                                end_q = int(rows[i]["Start_qry"])
                                lack_qry += qry_len - (end_q - args.ext)

                            len_align_r = ref_len - lack_ref
                            len_align_q = qry_len - lack_qry

                            identity = "/"
                            cov_r = "/"
                            cov_q = "/"

                            # Assign a quality score.
                            ## The assembled seq matches to the whole ref seq
                            if int(len_align_q) == ref_len:
                                quality_rq = 'A'
                            ## The assembled seq matches to the ref seq +-10% of ref length
                            elif int(len_align_q) in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
                                quality_rq = 'B'
                            ## The assembled seq matches to the ref seq, but not along all their length (>= 50% of their length align)
                            elif int(len_align_q) in range((ref_len - error_50_perc), (ref_len + error_50_perc)):
                                quality_rq = 'C'
                            else:
                                quality_rq = 'D'

                            # Write stats results in output file.
                            ## Local assembly performed with the DBG algorithm
                            if ".k" in qry_file.split('/')[-1]:
                                stats = [qry_id, g, flank, occ, k, a, strand, solution, len_q, ref, len_r, \
                                    start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                            ## Local assembly performed with the IRO algorithm
                            if ".dmax" in qry_file.split('/')[-1]:
                                stats = [qry_id, g, flank, occ, s, o, a, d, len_q, ref, len_r, \
                                    start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                            with open(ref_qry_sorted, "a") as output_ref:
                                output_ref.write('\n' + '\t'.join(str(i) for i in stats))

    except Exception as e:
        print("\nFile 'stats_alignment.py': Something wrong with the quality estimation for alignment in multiple chunks/flanks (DBG), when ref = reference sequence")
        print("Exception-")
        print(e)
        sys.exit(1)


#----------------------------------------------------
# Remove raw files
#----------------------------------------------------
# Remove the raw files obtained from statistics ('.log', '.delta', '.coords' files).
subprocess.run(["rm", nucmerLog])
subprocess.run(["rm", delta_file])
subprocess.run(["rm", coords_file])
subprocess.run(["rm", coords_sorted_file])

