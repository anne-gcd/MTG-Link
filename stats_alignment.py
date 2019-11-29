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
parser.add_argument("-ref", "--reference", action="store", help="file containing the reference sequence of the simulated gap (format: 'xxx.ingap.fasta') OR file containing the sequences of the scaffolds (format: 'xxx.fasta')", required=True)
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
    # Statistics about the alignment if Ref = reference sequence of simulated gap
    #-----------------------------------------------------------------------------
    if re.match('^.*.ingap.fasta$', args.reference):

        #Run NUCmer
        qry_id = qry_file.split('.k')[0]
        id_ = qry_id.split('gfa.')[1]
        prefix = id_ + ".ref_qry"

        nucmerLog = "{}_nucmer.log".format(id_)
        delta_file = prefix + ".delta"
        coords_file = prefix + ".coords"

        nucmer_command = ["nucmer", "-p", prefix, ref_file, qry_file]
        coords_command = ["show-coords", "-rcdlT", delta_file]

        with open(coords_file, "w") as coords, open(nucmerLog, "a") as log:
            subprocess.run(nucmer_command, stderr=log)
            subprocess.run(coords_command, stdout=coords, stderr=log)
            print(subprocess.check_output(coords_command))

        #Remove the raw file obtained from NUCmer
        if os.path.getsize(nucmerLog) <= 0:
            subprocess.run(["rm", nucmerLog])

        #Get scaffolds, kmer size, abundance min, gap size and chunk size values
        scaffolds = id_.split('.')[0]
        kmer_size = qry_file.split('.')[-5]
        k = int("".join(list(kmer_size)[1:]))
        abundance_min = qry_file.split('.')[-4]
        a = int("".join(list(abundance_min)[1:]))
        gap_size = id_.split('.')[1]
        g = int("".join(list(gap_size)[1:]))
        chunk_size = id_.split('.')[2]
        c = int("".join(list(chunk_size)[1:]))

        #Output stats file
        output_file = outDir + "/" + args.prefix + ".alignment.stats"
        stats_legend = ["Scaffolds", "Gap", "Chunk", "k", "a", "Strand", "Uniq_sol", "Len_R", "Len_Q", "Global_align", "Score", \
                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q"]

        #Alignment stats
        with open(ref_file, "r") as reference:
            for ref_record in SeqIO.parse(reference, "fasta"): #one loop (one reference seq)
                ref_seq = ref_record.seq
                ref_len = len(ref_seq)

        with open(qry_file, "r") as query:
            for qry_record in SeqIO.parse(query, "fasta"): #x records loops (x = nb of query (e.g. nb of inserted seq))
                qry_seq_id = qry_record.description.split(" ")[0]
                qry_seq = qry_record.seq
                qry_len = len(qry_seq) - k*2

                strand = "fwd"
                unique_sol = True
                global_align = "N"

                if "bkpt2" in str(qry_record.id):
                    strand = "rev"
                    ref_seq = rc(ref_seq)

                if "solution" in qry_record.description:
                    unique_sol = qry_record.description.split(" ")[-1]

                if qry_len > 0.9*ref_len and qry_len < 1.1*ref_len:
                    global_align = "Y"

                #Run Bio.Align
                alignments = aligner.align(qry_seq, ref_seq)
                score = alignments[0].score

                #Get output values from NUCmer
                reader = csv.DictReader(open(coords_file), \
                                        fieldnames=("S1", 'E1', "S2", "E2", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q", "COV_R", "COV_Q", "FRM_R", "FRM_Q", "TAG_1", "TAG_2"), \
                                        delimiter='\t')
                for row in reader:
                    if row["TAG_2"] == qry_seq_id:
                        len_r = row["LEN_R"]
                        len_q = row["LEN_Q"]
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

                        #Write stats results in output file
                        stats = [scaffolds, g, c, k, a, strand, unique_sol, len_r, len_q, global_align, score, \
                                start_r, end_r, start_q, end_q, len_align_r, len_align_q, identity, cov_r, cov_q, frame_r, frame_q]

                        if os.path.exists(output_file):
                            with open(output_file, "a") as output:
                                output.write('\n' + '\t'.join(str(i) for i in stats))
                        else:
                            with open(output_file, "a") as output:
                                output.write('\t'.join(j for j in stats_legend))
                                output.write('\n'+'\n' + '\t'.join(str(i) for i in stats))

                #Set back the reference sequence to the original one (not reversed)
                if strand == 'rev':
                    ref_seq = rc(ref_seq)
   

    #-----------------------------------------------------------------------------
    # Statistics about the alignment if Ref = scaffolds' sequences
    #-----------------------------------------------------------------------------
    if re.match('^.*.fasta$', args.reference):

        #----------------------------------------------------
        # Alignment Ref vs Qry (extension vs scaffolds' seq)
        #----------------------------------------------------
        #Run NUCmer to obtain alignment of extension portions (-ext) against the scaffolds' sequences
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


        #----------------------------------------------------
        # Alignment Qry vs Qry (fwd vs rev)
        #----------------------------------------------------
        #Run NUCmer to obtain alignment of fwd strand solution against rev strand solution
        prefix_qry = id_ + ".qry_qry"
        nucmerLog_qry = "{}_nucmer_qry_qry.log".format(id_)
        delta_file_qry = prefix_qry + ".delta"
        coords_file_qry = prefix_qry + ".coords"          

        nucmer_command_qry = ["nucmer", "-p", prefix_qry, qry_file, qry_file]
        coords_command_qry = ["show-coords", "-rcdlT", delta_file_qry]

        with open(coords_file_qry, "w") as coords_qry, open(nucmerLog_qry, "a") as log:
            subprocess.run(nucmer_command_qry, stderr=log)
            subprocess.run(coords_command_qry, stdout=coords_qry, stderr=log)

        #Output stats file of alignment query vs ref
        qry_qry_output = outDir + "/" + args.prefix + ".qry_qry.alignment.stats"
        stats_legend_qry = ["Gap", "Len_gap", "Chunk", "k", "a", "Solution1", "Len_Q1", "Solution2", "Len_Q2", \
                            "Start_Q1", "End_Q1", "Start_Q2", "End_Q2", "Len_align_Q1", "Len_align_Q2", "%_Id", "%_Cov_Q1", "%_Cov_Q2", "Frame_Q1", "Frame_Q2", "Quality"]

        #Get output values from NUCmer:
        reader_qry = csv.DictReader(open(coords_file_qry), \
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


except Exception as e:
    print("\nException-")
    print(e)
    sys.exit(1)