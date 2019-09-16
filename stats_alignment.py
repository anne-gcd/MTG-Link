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
parser = argparse.ArgumentParser(prog="stats_alignment.py", usage="%(prog)s -qry <query_sequences_file> -ref <reference_sequence> -p <output_file_prefix> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Statistics of the alignment of the inserted sequence obtained from MindTheGap (-qry) and the simulated gap (-ref)
                                Note: there are kmer flanking regions on the edges of the inserted sequence
                                '''))

parser.add_argument("-qry", "--query", action="store", help="file containing the inserted sequences obtained from MindTheGap (format: 'xxx.insertions.fasta')", required=True)
parser.add_argument("-ref", "--reference", action="store", help="file containing the reference sequence of the simulated gap (format: 'xxx.ingap.fasta')", required=True)
parser.add_argument("-p", "--prefix", action="store", help="prefix of output file to save the statistical results", required=True)
parser.add_argument("-out", "--outDir", action="store", default="./mtg10x_results/alignments_stats", help="output directory for saving results")

args = parser.parse_args()

if re.match('^.*.insertions.fasta$', args.query) is None:
    parser.error("The suffix of the inserted sequences (query sequences) file should be: '.insertions.fasta'")

if re.match('^.*.ingap.fasta$', args.reference) is None:
    parser.error("The suffix of the reference sequence file should be: '.ingap.fasta'")

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

#----------------------------------------------------
# Statistics about the alignment
#----------------------------------------------------
try:
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
        for qry_record in SeqIO.parse(query, "fasta"): #x records loops (x = nb of query/inserted seq)
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
   

except Exception as e:
    print("\nException-")
    print(e)
    sys.exit(1)