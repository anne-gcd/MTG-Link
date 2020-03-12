#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import re
import csv
import argparse
import subprocess
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO, Align
from Bio.Seq import Seq


#PairwiseAligner object
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
                                Statistics about the inserted sequence obtained from MindTheGap (-qry)
                                Note: there are kmer flanking regions on the edges of the inserted sequence
                                '''))

parser.add_argument("-qry", "--query", action="store", help="file containing the inserted sequences obtained from MindTheGap (format: 'xxx.insertions.fasta')", required=True)
parser.add_argument("-ref", "--reference", action="store", help="file containing the reference sequence of the simulated gap (format: 'xxx.ingap.fasta') or the sequences of the contigs (format: 'xxx.fasta')", required=True)
parser.add_argument("-ext", "--ext", action="store", type=int, help="size of the gap, on both sides; determine start/end of gapfilling")
parser.add_argument("-p", "--prefix", action="store", help="prefix of output file to save the statistical results", required=True)
parser.add_argument("-out", "--outDir", action="store", default="./mtg10x_results/alignments_stats", help="output directory for saving results")

args = parser.parse_args()

if re.match('^.*.insertions.fasta$', args.query) is None:
    parser.error("The suffix of the inserted sequences (query sequences) file should be: '.insertions.fasta'")

if re.match('^.*.fasta$', args.reference) is None:
    parser.error("The suffix of the reference sequence file should be: '.fasta'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
qry_file = os.path.abspath(args.query)
if not os.path.exists(args.query):
    parser.error("The path of the query file (inserted sequences file) doesn't exist")
print("\nQuery file: " + qry_file)

ref_file = os.path.abspath(args.reference)
if not os.path.exists(ref_file):
    parser.error("The path of the reference file doesn't exist")
print("Reference file: " + ref_file)

#----------------------------------------------------
# Directory for saving results
#----------------------------------------------------
cwd = os.getcwd()
if not os.path.exists(args.outDir):
    os.mkdir(args.outDir)
try:
    os.chdir(args.outDir)
except:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()
print("\nThe results are saved in " + outDir)

try:
    #-----------------------------------------------------------------------------
    # Statistics about the Alignment Ref vs Qry
    #-----------------------------------------------------------------------------

    if re.match('^.*.ingap.fasta$', args.reference):
        #----------------------------------------------------
        # Ref = reference sequence of simulated gap
        #----------------------------------------------------
        #Run NUCmer to obtain alignment of extension portions (-ext) against the contigs' sequences
        qry_id = qry_file.split('.')[-9]
        qry_k = qry_file.split('.')[-6]
        qry_a = qry_file.split('.')[-5]
        id_ = str(qry_id) + "." + str(qry_k) + "." + str(qry_a)
        prefix = id_ + ".ref_qry"

        nucmerLog = "{}_nucmer_ref_qry.log".format(id_)
        delta_file = prefix + ".delta"
        coords_file = prefix + ".coords.unsorted" 

        nucmer_command = ["nucmer", "--maxmatch", "-p", prefix, ref_file, qry_file]
        coords_command = ["show-coords", "-rcdlT", delta_file]

        with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
            subprocess.run(nucmer_command, stderr=log)
            subprocess.run(coords_command, stdout=coords, stderr=log)
        
        #Sort the 'xxx.coords.unsorted' file for further analysis
        coords_sorted_file = prefix + ".coords"
        sort_command = ["sort", "-n", coords_file]
        with open(coords_sorted_file, "w") as coords_sorted:
            subprocess.run(sort_command, stdout=coords_sorted)

        #Output stats file of alignment query vs ref
        ref_qry_output = outDir + "/" + args.prefix + ".ref_qry.alignment.stats.unsorted"
        stats_legend = ["Gap", "Len_gap", "Chunk", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]
        
        #Get the gap size, chunk size, kmer size and abundance min values
        gap_size = qry_file.split('.')[-8]
        g = int("".join(list(gap_size)[1:]))
        chunk_size = qry_file.split('.')[-7]
        c = int("".join(list(chunk_size)[1:]))
        k = int("".join(list(qry_k)[1:]))
        a = int("".join(list(qry_a)[1:]))

        #Get output values from NUCmer:
        reader = csv.DictReader(open(coords_sorted_file), \
                                    fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                    delimiter='\t')

        rows = list(reader)
        for row in rows[3:]:

            len_q = row["LEN_Q"]
            ref = row["TAG_1"]
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
            solution = str(row["TAG_2"]).split('_sol_')[1]

            if "bkpt1" in str(row["TAG_2"]):
                strand = "fwd"
            else:
                strand = "rev"

            # Estimate quality of gapfilled sequence
            ref_len = int(len_r)
            error_10_perc = int(0.1 * ref_len)
            qry_len = int(len_q) - 2*args.ext

            #lengths of both sequences are equal +-10% of ref length
            if qry_len in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
                #the gapfilled seq matches to the whole ref seq
                if int(len_align_q) == ref_len:
                    quality_rq = 'A'
                #the gapfilled seq matches to the ref seq +-10% of ref length
                elif int(len_align_q) in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
                    quality_rq = 'B'
                #the gapfilled seq matches to the ref seq, but not along all their length (>= 50% of their length align)
                elif int(len_align_q) >= int(0.5*ref_len) and int(len_align_r) >= int(0.5*qry_len):
                    quality_rq = 'C'
                else:
                    quality_rq = 'D'

            else:
                quality_rq = 'D'

            #Write stats results in output file
            stats = [qry_id, g, c, k, a, strand, solution, len_q, ref, len_r, \
                    start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

            if os.path.exists(ref_qry_output):
                with open(ref_qry_output, "a") as output:
                    output.write('\n' + '\t'.join(str(i) for i in stats))
            else:
                with open(ref_qry_output, "a") as output:
                    output.write('\t'.join(j for j in stats_legend))
                    output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))

        #Sort the 'xxx.alignment.stats.unsorted' file for further analysis
        ref_qry_sorted = outDir + "/" + args.prefix + ".ref_qry.alignment.stats"
        order_command = ["sort", "-k6,7", "-k11,12n", "-r", ref_qry_output]
        with open(ref_qry_sorted, "w") as r_sorted:
            subprocess.run(order_command, stdout=r_sorted)

        #If alignment in multiple chunks, calculate the appropriate quality score
        with open(ref_qry_sorted, "r") as r:
            r.seek(0)
            reader = csv.DictReader(r, fieldnames=("Gap", "Len_gap", "Chunk", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                                    "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
                                    delimiter='\t')

            rows = list(reader)
            for i in range(1, len(rows)):
                if ("fwd" in rows[i]["Strand"]) or ("rev" in rows[i]["Strand"]):

                    qry_id = rows[i]["Gap"]
                    g = rows[i]["Len_gap"]
                    c = rows[i]["Chunk"]
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

                    if (strand != rows[i-1]["Strand"]) or (solution != rows[i-1]["Solution"]):
                        lack_ref = 0
                        lack_qry = 0
                        start_r = int(rows[i]["Start_ref"])
                        lack_ref += start_r - 1
                        if frame_q == "1":
                            start_q = int(rows[i]["Start_qry"])
                            lack_qry += start_q - (args.ext+1)
                        elif frame_q == "-1":
                            end_q = int(rows[i]["Start_qry"])
                            lack_qry += qry_len - (end_q - args.ext)

                    if (strand == rows[i-1]["Strand"]) and (solution == rows[i-1]["Solution"]):
                        if int(rows[i]["Start_ref"]) > int(rows[i-1]["End_ref"]):
                            lack_ref += int(rows[i]["Start_ref"]) - int(rows[i-1]["End_ref"])
                        if frame_q == "1":
                            if int(rows[i]["Start_qry"]) > int(rows[i-1]["End_qry"]):
                                lack_qry += int(rows[i]["Start_qry"]) - int(rows[i-1]["End_qry"])
                        elif frame_q == "-1":
                            if int(rows[i]["Start_qry"]) < int(rows[i-1]["End_qry"]):
                                lack_qry += int(rows[i-1]["End_qry"]) - int(rows[i]["Start_qry"])

                    if (i == len(rows)-1) or ((strand != rows[i+1]["Strand"]) or (solution != rows[i+1]["Solution"])):
                        end_r = int(rows[i]["End_ref"])
                        lack_ref += ref_len - end_r
                        if frame_q == "1":
                            end_q = int(rows[i]["End_qry"])
                            lack_qry += qry_len - (end_q - args.ext)
                        elif frame_q == "-1":
                            start_q = int(rows[i]["End_qry"])
                            lack_qry += start_q - (args.ext+1)

                        len_align_r = ref_len - lack_ref
                        len_align_q = qry_len - lack_qry

                        identity = "/"
                        cov_r = "/"
                        cov_q = "/"

                        #Assign a quality score
                        #the gapfilled seq matches to the whole ref seq
                        if int(len_align_q) == ref_len:
                            quality_rq = 'A'
                        #the gapfilled seq matches to the ref seq +-10% of ref length
                        elif int(len_align_q) in range((ref_len - error_10_perc), (ref_len + error_10_perc)):
                            quality_rq = 'B'
                        #the gapfilled seq matches to the ref seq, but not along all their length (>= 50% of their length align)
                        elif int(len_align_q) >= int(0.5*ref_len) and int(len_align_r) >= int(0.5*qry_len):
                            quality_rq = 'C'
                        else:
                            quality_rq = 'D'

                        #Write stats results in output file
                        stats = [qry_id, g, c, k, a, strand, solution, len_q, ref, len_r, \
                                start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                        with open(ref_qry_sorted, "a") as output_ref:
                            output_ref.write('\n' + '\t'.join(str(i) for i in stats))


    if re.match('^.*.fasta$', args.reference):
        #----------------------------------------------------
        # Ref = contigs' sequences
        #----------------------------------------------------
        #Run NUCmer to obtain alignment of extension portions (-ext) against the contigs' sequences
        qry_id = qry_file.split('.')[-9]
        qry_k = qry_file.split('.')[-6]
        qry_a = qry_file.split('.')[-5]
        id_ = str(qry_id) + "." + str(qry_k) + "." + str(qry_a)
        prefix = id_ + ".ref_qry"

        nucmerLog = "{}_nucmer_ref_qry.log".format(id_)
        delta_file = prefix + ".delta"
        coords_file = prefix + ".coords"

        nucmer_command = ["nucmer", "-p", prefix, ref_file, qry_file]
        coords_command = ["show-coords", "-rcdlT", delta_file]

        with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
            subprocess.run(nucmer_command, stderr=log)
            subprocess.run(coords_command, stdout=coords, stderr=log)

        #Output stats file of alignment query vs ref
        ref_qry_output = outDir + "/" + args.prefix + ".ref_qry.alignment.stats"
        stats_legend = ["Gap", "Len_gap", "Chunk", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"]
        
        #Get the gap size, chunk size, kmer size and abundance min values
        gap_size = qry_file.split('.')[-8]
        g = int("".join(list(gap_size)[1:]))
        if g == 0:
            g = "NA"
        chunk_size = qry_file.split('.')[-7]
        c = int("".join(list(chunk_size)[1:]))
        k = int("".join(list(qry_k)[1:]))
        a = int("".join(list(qry_a)[1:]))

        #Get output values from NUCmer:
        reader = csv.DictReader(open(coords_file), \
                                    fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                    delimiter='\t')

        rows = list(reader)
        for row in rows[3:]:
            if row["TAG_1"] in str(qry_id):

                len_q = row["LEN_Q"]
                ref = row["TAG_1"]
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
                solution = str(row["TAG_2"]).split('_sol_')[1]

                if "bkpt1" in str(row["TAG_2"]):
                    strand = "fwd"
                else:
                    strand = "rev"

                # Estimate quality of gapfilled sequence
                left = str(qry_id).split('_')[0]
                left_scaffold = left[:-1]
                right = str(qry_id).split('_')[1]
                right_scaffold = right[:-1]
                error_10_perc = int(0.1 * args.ext)

                #ref = Left scaffold
                if ref == left_scaffold:
                    #extension of qry match perfectly as expected to ref
                    if (('+' in left and int(start_r) == (int(len_r) - args.ext + 1) and int(end_r) == int(len_r)) and ((strand == "fwd" and int(start_q) == 1 and int(end_q) == args.ext) or (strand == "rev" and int(start_q) == int(len_q) and int(end_q) == (int(len_q) - args.ext + 1)))) \
                        or (('-' in left and int(start_r) == 1 and int(end_r) == args.ext) and ((strand == "fwd" and int(start_q) == args.ext and int(end_q) == 1) or (strand == "rev" and int(start_q) == (int(len_q) - args.ext + 1) and int(end_q) == int(len_q)))):
                        quality_rq = 'A'
                    #extension of qry almost match as expected to ref (+-10% of extension size) 
                    elif (('+' in left and int(start_r) in range((int(len_r)-args.ext+1 - error_10_perc), (int(len_r)-args.ext+1 + error_10_perc+1)) and int(end_r) in range((int(len_r) - error_10_perc), (int(len_r)+1))) and ((strand == "fwd" and int(start_q) in range(1, (1 + error_10_perc+1)) and int(end_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1))) or (strand == "rev" and int(start_q) in range((int(len_q) - error_10_perc), (int(len_q)+1)) and int(end_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1))))) \
                        or (('-' in left and int(start_r) in range(1, (1 + error_10_perc+1)) and int(end_r) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1))) and ((strand == "fwd" and int(start_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1)) and int(end_q) in range(1, (1 + error_10_perc+1))) or (strand == "rev" and int(start_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1)) and int(end_q) in range((int(len_q) - error_10_perc), (int(len_q)+1))))):
                        quality_rq = 'B'
                    #extension of qry almost match (+-ext) as expected to ref
                    elif (('+' in left and int(start_r) >= (int(len_r) - args.ext + 1)) and ((strand == "fwd" and int(end_q) <= 2*args.ext) or (strand == "rev" and int(end_q) >= (int(len_q) - 2*args.ext)))) \
                        or (('-' in left and int(end_r) <= args.ext) and ((strand == "fwd" and int(start_q) <= 2*args.ext) or (strand == "rev" and int(start_q) >= (int(len_q) - 2*args.ext)))):
                        quality_rq = 'C'
                    else:
                        quality_rq = 'D'

                #ref = Right scaffold
                elif ref == right_scaffold:
                    #extension of qry match perfectly as expected to ref
                    if (('+' in right and int(start_r) == 1 and int(end_r) == args.ext) and ((strand == "fwd" and int(start_q) == (int(len_q) - args.ext + 1) and int(end_q) == int(len_q)) or (strand == "rev" and int(start_q) == args.ext and int(end_q) == 1))) \
                        or (('-' in right and int(start_r) == (int(len_r) - args.ext + 1) and int(end_r) == int(len_r)) and ((strand == "fwd" and int(start_q) == int(len_q) and int(end_q) == (int(len_q) - args.ext +1)) or (strand == "rev" and int(start_q) == 1 and int(end_q) == args.ext))):
                        quality_rq = 'A'
                    #extension of qry almost match as expected to ref (+-10% of extension size) 
                    elif (('+' in right and int(start_r) in range(1, (1 + error_10_perc+1)) and int(end_r) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1))) and ((strand == "fwd" and int(start_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1)) and int(end_q) in range((int(len_q) - error_10_perc), (int(len_q)+1))) or (strand == "rev" and int(start_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1)) and int(end_q) in range(1, (1 + error_10_perc+1))))) \
                        or (('-' in right and int(start_r) in range((int(len_r)-args.ext+1 - error_10_perc), (int(len_r)-args.ext+1 + error_10_perc+1)) and int(end_r) in range((int(len_r) - error_10_perc), (int(len_r)+1))) and ((strand == "fwd" and int(start_q) in range((int(len_q) - error_10_perc), (int(len_q)+1)) and int(end_q) in range((int(len_q)-args.ext+1 - error_10_perc), (int(len_q)-args.ext+1 + error_10_perc+1))) or (strand == "rev" and int(start_q) in range(1, (1 + error_10_perc+1)) and int(end_q) in range((args.ext - error_10_perc), (args.ext + error_10_perc+1))))):
                        quality_rq = 'B'
                    #extension of qry almost match (+-ext) as expected to ref
                    elif (('+' in right and int(end_r) <= args.ext) and ((strand == "fwd" and int(start_q) >= (int(len_q) - 2*args.ext)) or (strand == "rev" and int(start_q) <= 2*args.ext))) \
                        or (('-' in right and int(start_r) >= (int(len_r) - args.ext + 1)) and ((strand == "fwd" and int(end_q) >= (int(len_q) - 2*args.ext)) or (strand == "rev" and int(end_q) <= 2*args.ext))):
                        quality_rq = 'C'
                    else:
                        quality_rq = 'D'
                
                #Write stats results in output file
                stats = [qry_id, g, c, k, a, strand, solution, len_q, ref, len_r, \
                        start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q, quality_rq]

                if os.path.exists(ref_qry_output):
                    with open(ref_qry_output, "a") as output:
                        output.write('\n' + '\t'.join(str(i) for i in stats))
                else:
                    with open(ref_qry_output, "a") as output:
                        output.write('\t'.join(j for j in stats_legend))
                        output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))


    #-----------------------------------------------------------------------------
    # Statistics about the Alignment Qry vs Qry (fwd vs rev)
    #-----------------------------------------------------------------------------
    #Run NUCmer to obtain alignment of fwd strand solution against rev strand solution
    prefix_qry = id_ + ".qry_qry"
    nucmerLog_qry = "{}_nucmer_qry_qry.log".format(id_)
    delta_file_qry = prefix_qry + ".delta"
    coords_file_qry = prefix_qry + ".coords.unsorted"          

    nucmer_command_qry = ["nucmer", "--maxmatch", "-r", "-p", prefix_qry, qry_file, qry_file]
    coords_command_qry = ["show-coords", "-rcdlT", delta_file_qry]

    with open(coords_file_qry, "w") as coords_qry, open(nucmerLog_qry, "a") as log:
        subprocess.run(nucmer_command_qry, stderr=log)
        subprocess.run(coords_command_qry, stdout=coords_qry, stderr=log)
    
    #Sort the 'xxx.coords.unsorted' file for further analysis
    coords_qry_sorted_file = prefix_qry + ".coords"
    sort_command_qry = ["sort", "-n", coords_file_qry]
    with open(coords_qry_sorted_file, "w") as coords_qry_sorted:
        subprocess.run(sort_command_qry, stdout=coords_qry_sorted)

    #Output stats file of alignment query vs ref
    qry_qry_output = outDir + "/" + args.prefix + ".qry_qry.alignment.stats.unsorted"
    stats_legend_qry = ["Gap", "Len_gap", "Chunk", "k", "a", "Solution1", "Len_Q1", "Solution2", "Len_Q2", \
                        "Start_Q1", "End_Q1", "Start_Q2", "End_Q2", "Len_align_Q1", "Len_align_Q2", "%_Id", "%_Cov_Q1", "%_Cov_Q2", "Frame_Q1", "Frame_Q2", "Quality"]

    #Get output values from NUCmer:
    reader_qry = csv.DictReader(open(coords_qry_sorted_file), \
                                fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                delimiter='\t')

    rows_qry = list(reader_qry)
    for row in rows_qry[3:]:

        len_q1 = row["LEN_R"]
        len_q2 = row["LEN_Q"]
        start_q1 = row["S1"]
        end_q1 = row["E1"]
        start_q2 = row["S2"]
        end_q2 = row["E2"]
        len_align_q1 = row["LEN_1"]
        len_align_q2 = row["LEN_2"]
        identity = row["%_IDY"]
        cov_q1 = row["COV_R"]
        cov_q2 = row["COV_Q"]
        frame_q1 = row["FRM_R"]
        frame_q2 = row["FRM_Q"]

        if "bkpt1" in str(row["TAG_1"]):
            strand_q1 = "fwd"
        else:
            strand_q1 = "rev"

        if "bkpt1" in str(row["TAG_2"]):
            strand_q2 = "fwd"
        else:
            strand_q2 = "rev"
        
        solution_1 = strand_q1 + str(row["TAG_1"]).split('_sol_')[1]
        solution_2 = strand_q2 + str(row["TAG_2"]).split('_sol_')[1]

        # Estimate quality of gapfilled sequence (only for fwd vs rc_rev)
        if solution_1 != solution_2:
            
            #lengths of both sequences are equal +-10%
            if int(len_q1) in range((int(len_q2) - int(0.1*int(len_q2))), (int(len_q2) + int(0.1*int(len_q2)))):
                #both sequences match to each other along all their length
                if int(len_align_q1) == int(len_q1) and int(len_align_q2) == int(len_q2):
                    quality_qq = 'A'
                #one sequence match to the other along all its length or both sequences match to each other along their length +-10%
                elif (int(len_align_q1) == int(len_q1) or int(len_align_q2) == int(len_q2)) or (int(len_align_q1) in range((int(len_q1) - int(0.1*int(len_q1))), (int(len_q1) + int(0.1*int(len_q1)))) or (int(len_align_q2) in range((int(len_q2) - int(0.1*int(len_q2))), (int(len_q2) + int(0.1*int(len_q2)))))):
                    quality_qq = 'B'
                #both sequences match to each other, but not along all their length (>= 50% of their length align)
                elif int(len_align_q1) >= int(0.5*int(len_q1)) and int(len_align_q2) >= int(0.5*int(len_q2)):
                    quality_qq = 'C'
                else:
                    quality_qq = 'D'
            #lengths of both sequences are not equal +-10%
            else:
                quality_qq = 'D'

        else:
            quality_qq = 'D'

        #Write stats results in output file
        stats_qry = [qry_id, g, c, k, a, solution_1, len_q1, solution_2, len_q2, \
                    start_q1, end_q1, start_q2, end_q2, len_align_q1, len_align_q2, identity, cov_q1, cov_q2, frame_q1, frame_q2, quality_qq]

        if os.path.exists(qry_qry_output):
            with open(qry_qry_output, "a") as output_qry:
                output_qry.write('\n' + '\t'.join(str(i) for i in stats_qry))
        else:
            with open(qry_qry_output, "a") as output_qry:
                output_qry.write('\t'.join(j for j in stats_legend_qry))
                output_qry.write('\n'+'\n' + '\t'.join(str(i) for i in stats_qry))

    #sort the 'xxx.align;ent.stats.unsorted' file for further analysis
    qry_qry_sorted = outDir + "/" + args.prefix + ".qry_qry.alignment.stats"
    order_command_qry = ["sort", "-n", qry_qry_output]
    with open(qry_qry_sorted, "w") as f_sorted:
        subprocess.run(order_command_qry, stdout=f_sorted)

    #If alignment in multiple chunks, calculate the appropriate quality score
    with open(qry_qry_sorted, "r") as f:
        f.seek(0)
        reader = csv.DictReader(f, fieldnames=("Gap", "Len_gap", "Chunk", "k", "a", "Solution1", "Len_Q1", "Solution2", "Len_Q2", \
                                "Start_Q1", "End_Q1", "Start_Q2", "End_Q2", "Len_align_Q1", "Len_align_Q2", "%_Id", "%_Cov_Q1", "%_Cov_Q2", "Frame_Q1", "Frame_Q2", "Quality"), \
                                delimiter='\t')

        rows = list(reader)
        lack_Q1 = 0
        lack_Q2 = 0
        for i in range(1, len(rows_qry)):
            if ("fwd" in rows[i]["Solution1"]) and ("rev" in rows[i]["Solution2"]):

                qry_id = rows[i]["Gap"]
                g = rows[i]["Len_gap"]
                c = rows[i]["Chunk"]
                k = rows[i]["k"]
                a = rows[i]["a"]
                solution1 = rows[i]["Solution1"]
                solution2 = rows[i]["Solution2"]
                len_Q1 = int(rows[i]["Len_Q1"])
                len_Q2 = int(rows[i]["Len_Q2"])
                frame_Q1 = rows[i]["Frame_Q1"]
                frame_Q2 = rows[i]["Frame_Q2"]

                if (solution1 != rows[i-1]["Solution1"]) or (solution2 != rows[i-1]["Solution2"]):
                    start_Q1 = int(rows[i]["Start_Q1"])
                    end_Q2 = int(rows[i]["Start_Q2"])
                    lack_Q1 += (start_Q1 - 1)
                    lack_Q2 += (len_Q2 - end_Q2)

                if (solution1 == rows[i-1]["Solution1"]) and (solution2 == rows[i-1]["Solution2"]):
                    if int(rows[i]["Start_Q1"]) > int(rows[i-1]["End_Q1"]):
                        lack_Q1 += (int(rows[i]["Start_Q1"]) - int(rows[i-1]["End_Q1"]))

                    if int(rows[i]["Start_Q2"]) < int(rows[i-1]["End_Q2"]):
                        lack_Q2 += (int(rows[i-1]["End_Q2"]) - int(rows[i]["Start_Q2"]))

                if (solution1 != rows[i+1]["Solution1"]) or (solution2 != rows[i+1]["Solution2"]):
                    end_Q1 = int(rows[i]["End_Q1"])
                    start_Q2 = int(rows[i]["End_Q2"])
                    lack_Q1 += (len_Q1 - end_Q1)
                    lack_Q2 += (start_Q2 - 1)

            elif ("fwd" in rows[i-1]["Solution1"]) and ("rev" in rows[i-1]["Solution2"]):
                break

        len_align_Q1 = len_Q1 - lack_Q1
        len_align_Q2 = len_Q2 - lack_Q2

        identity = "/"
        cov_Q1 = "/"
        cov_Q2 = "/"

        #Assign a quality score
        #both sequences match to each other along all their length
        if len_align_Q1 == len_Q1 and len_align_Q2 == len_Q2:
            quality_qq = 'A'
        #one sequence match to the other along all its length or both sequences match to each other along their length +-10%
        elif (len_align_Q1 == len_Q1 or len_align_Q2 == len_Q2) or (len_align_Q1 in range((len_Q1 - int(0.1*len_Q1)), (len_Q1 + int(0.1*len_Q1))) or len_align_Q2 in range((len_Q2 - int(0.1*len_Q2)), (len_Q2 + int(0.1*len_Q2)))):
            quality_qq = 'B'
        #both sequences match to each other, but not along all their length (>= 50% of their length align)
        elif len_align_Q1 >= int(0.5*len_Q1) and len_align_Q2 >= int(0.5*len_Q2):
            quality_qq = 'C'
        else:
            quality_qq = 'D'

        #Write stats results in output file
        stats_qry = [qry_id, g, c, k, a, solution1, len_Q1, solution2, len_Q2, \
                    start_Q1, end_Q1, start_Q2, end_Q2, len_align_Q1, len_align_Q2, identity, cov_Q1, cov_Q2, frame_Q1, frame_Q2, quality_qq]

        with open(qry_qry_sorted, "a") as output_qry:
            output_qry.write('\n' + '\t'.join(str(i) for i in stats_qry))


    #Remove the raw file obtained from statistics ('.delta', '.coords', '.unsorted' files)
    subprocess.run(["rm", delta_file])
    subprocess.run(["rm", delta_file_qry])
    subprocess.run(["rm", coords_file])
    subprocess.run(["rm", coords_file_qry])


except Exception as e:
    print("\nException-")
    print(e)
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    sys.exit(1)