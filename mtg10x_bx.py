#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import argparse
import csv
import re
import subprocess
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO, Align
from helpers_bx import Gap, Scaffold, extract_barcodes, get_reads, mtg_fill, stats_align, get_position_for_edges, output_gfa_with_solution


#----------------------------------------------------
# Arg parser
#----------------------------------------------------
parser = argparse.ArgumentParser(prog="mtg10x.py", usage="%(prog)s -in <GFA_file> -c <chunk_size> -bam <BAM_file> -reads <reads_file> -index <index_file> -f <freq_barcodes> [options]", \
                                formatter_class=argparse.RawTextHelpFormatter, \
                                description=(''' \
                                Gapfilling with 10X data, using MindTheGap in 'breakpoint' mode
                                ----------------------------------------------------------------
                                To be able to execute this script, you need to install a virtual environment containing Samtools, Biopython and Gfapy.
                                You also need to install Con10X (especially BamExtractor and GetReads) and MindTheGap 
                                Use MTG in breakpoint mode, but taking an offset of size k
                                
                                [Main options]: '-gfa': input GFA file containing the contigs paths
                                                '-c': chunk size 
                                                '-bam': BAM file
                                                '-reads': file of indexed reads
                                                '-index': barcodes index file
                                                '-f': minimal frequence of extracted barcodes from BAM file [default '2']
                                                '-out': output directory [default './mtg10x_results']
                                                '-refDir': directory containing the reference sequences if any
                                                '-scaffolds': file containing the sequences of the scaffolds
                                                '-line': line of GFA file input from which to start analysis 

                                [MindTheGap options]: '-bkpt': breakpoint file (with possibly offset of size k removed)              
                                                      '-k': size of a kmer [default '[51, 41, 31, 21]']
                                                      '--force': to force search on all k values
                                                      '-a': minimal abundance threshold for solid kmers [default '[3, 2]']
                                                      '-ext': extension size of the gap on both sides
                                                      '-max-nodes': maximum number of nodes in contig graph [default '1000']
                                                      '-max-length': maximum length of gapfilling (nt) [default '10000']
                                                      '-nb-cores': number of cores [default '4']
                                                      '-max-memory': max memory for graph building (in MBytes) [default '8000']
                                                      '-verbose': verbosity level  [default '1']
                                '''))

parserMain = parser.add_argument_group("[Main options]")
parserMtg = parser.add_argument_group("[MindTheGap option]")

parserMain.add_argument('-in', action="store", dest="input", help="input GFA file containing the contigs paths (format: xxx.gfa)", required=True)
parserMain.add_argument('-c', action="store", dest="chunk", type=int, help="size of the chunk for gapfilling", required=True)
parserMain.add_argument('-bam', action="store", dest="bam", help="BAM file containing the reads of the individual (format: xxx.bam)", required=True)
parserMain.add_argument('-reads', action="store", dest="reads", help="file of indexed reads", required=True)
parserMain.add_argument('-index', action="store", dest="index", help="barcodes index file", required=True)
parserMain.add_argument('-f', action="store", dest="freq", type=int, default=2, help="minimal frequence of extracted barcodes from BAM file, in a specific region (chunk)")
parserMain.add_argument('-out', action="store", dest="out_dir", default="./mtg10x_results", help="output directory for result files")
parserMain.add_argument('-refDir', action="store", dest="refDir", help="directory containing the reference sequences if any, for statistical purposes")
parserMain.add_argument('-scaffolds', action="store", dest="scaffs", help="file containing the sequences of the scaffolds (format: xxx.fasta)", required=True)
parserMain.add_argument('-line', action="store", dest="line", type=int, help="line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input)")

parserMtg.add_argument('-bkpt', action="store", dest="bkpt", help="breakpoint file in fasta format")
parserMtg.add_argument('-k', action="store", dest= "k_mtg", default=[51, 41, 31, 21],  nargs='*', type=int, help="kmer size used for gapfilling")
parserMtg.add_argument("--force", action="store_true", help="to force search on all k values")
parserMtg.add_argument('-a', action="store", dest="a_mtg", default=[3, 2], nargs='*', type=int, help="minimal abundance of kmers used for gapfilling")
parserMtg.add_argument('-ext', action="store", dest="ext", type=int, help="size of the extension of the gap, on both sides [by default k]; determine start/end of gapfilling")
parserMtg.add_argument('-max-nodes', action="store", dest="max_nodes_mtg", type=int, default=1000, help="maximum number of nodes in contig graph")
parserMtg.add_argument('-max-length', action="store", dest="max_length_mtg", type=int, default=10000, help="maximum length of gapfilling (nt)")
parserMtg.add_argument('-nb-cores', action="store", dest="nb_cores_mtg", type=int, default=4, help="number of cores")
parserMtg.add_argument('-max-memory', action="store", dest="max_memory_mtg", type=int, default=8000, help="max memory for graph building (in MBytes)")
parserMtg.add_argument('-verbose', action="store", dest="verbose_mtg", type=int, default=1, help="verbosity level")

args = parser.parse_args()

if re.match('^.*.gfa$', args.input) is None:
    parser.error("The suffix of the GFA file should be: '.gfa'")

if re.match('^.*.bam$', args.bam) is None:
    parser.error("The suffix of the BAM file should be: '.bam'")

if re.match('^.*.fasta$', args.scaffs) is None:
    parser.error("The suffix of the file containing the sequences of the scaffolds should be: '.fasta'")

#----------------------------------------------------
# Input files
#----------------------------------------------------
gfa_file = os.path.abspath(args.input)
if not os.path.exists(gfa_file):
    parser.error("The path of the GFA file doesn't exist")
gfa_name = gfa_file.split('/')[-1]
print("\nInput GFA file: " + gfa_file)

bam_file = os.path.abspath(args.bam)
if not os.path.exists(bam_file): 
    parser.error("The path of the BAM file doesn't exist")
print("BAM file: " + bam_file)

reads_file = os.path.abspath(args.reads)
if not os.path.exists(reads_file):
    parser.error("The path of the file of indexed reads doesn't exist")
print("File of indexed reads: " + reads_file)

index_file = os.path.abspath(args.index)
print("Barcodes index file (prefix): " + index_file)

if args.refDir is not None:
    refDir = os.path.abspath(args.refDir)
    if not os.path.exists(refDir):
        parser.error("The path of the directory containing the reference sequences doesn't exist")

scaffs_file = os.path.abspath(args.scaffs)
if not os.path.exists(scaffs_file):
    parser.error("The path of the file of scaffolds' sequences doesn't exist")
print("File of scaffolds' sequences: " + scaffs_file)

#----------------------------------------------------
# Directories for saving results
#----------------------------------------------------
cwd = os.getcwd() 

#outDir
if not os.path.exists(args.out_dir):
    os.mkdir(args.out_dir)
try:
    os.chdir(args.out_dir)
except:
    print("Something wrong with specified directory. Exception-", sys.exc_info())
    print("Restoring the path")
    os.chdir(cwd)
outDir = os.getcwd()
print("\nThe results are saved in " + outDir)

#statsDir
statsDir = outDir + "/alignments_stats"

#----------------------------------------------------
# Gapfilling pipeline
#----------------------------------------------------
try:
    gfa = gfapy.Gfa.from_file(gfa_file)
    out_gfa_file = "mtg10x_" + gfa_name

    #----------------------------------------------------
    # GFA output: case no gap
    #----------------------------------------------------
    #If no gap, rewrite all the lines into GFA output
    if len(gfa.gaps) == 0:
        with open(out_gfa_file, "w") as f:
            out_gfa = gfapy.Gfa()
            for line in gfa.lines:
                out_gfa.add_line(str(line))
            out_gfa.to_file(out_gfa_file)

    #----------------------------------------------------   
    # Fill the gaps
    #----------------------------------------------------
    #If gap, rewrite the H and S lines into GFA output
    if args.line is None:
        with open(out_gfa_file, "w") as f:
            out_gfa = gfapy.Gfa()
            out_gfa.add_line("H\tVN:Z:2.0")
            for line in gfa.segments:
                out_gfa.add_line(str(line))
            out_gfa.to_file(out_gfa_file)
        
    gaps = []
    #If '-line' argument provided, start analysis from this line in GFA file input
    if args.line is not None:
        for _gap_ in gfa.gaps[(args.line - (len(gfa.segments)+2)):]:
            gaps.append(_gap_)
    else:
        for _gap_ in gfa.gaps:
            gaps.append(_gap_)

    for current_gap in gaps:
        gap = Gap(current_gap)
        gap.info()
        gap_label = gap.label()
        
        left_scaffold = Scaffold(current_gap, gap.left)
        right_scaffold = Scaffold(current_gap, gap.right)

        #If chunk size larger than length of scaffold(s), abort tentative of gapfilling for this scaffold
        if args.chunk > left_scaffold.len or args.chunk > right_scaffold.len:
            args.chunk = min(left_scaffold.len, right_scaffold.len)
            print("The chunk size was not smaller than the length of the scaffolds: it was set to the minimal scaffold length\n")

        #Save current G line into a temporary file
        with open("tmp.gap", "w") as tmp_gap:
            tmp_gap.write(str(current_gap))
            tmp_gap.seek(0)

        #----------------------------------------------------
        # BamExtractor
        #----------------------------------------------------
        #Initiate a dictionary to count the occurences of each barcode
        barcodes_occ = {}
        
        #Obtain the left barcodes and store the elements in a set
        left_region = left_scaffold.chunk(args.chunk)
        left_barcodes_file = "{}{}.c{}.left.barcodes".format(left_scaffold.name, left_scaffold.orient, args.chunk)

        print("Barcodes from left chunk ({}): {}".format(left_region, left_barcodes_file))
        with open(left_barcodes_file, "w+") as left_barcodes:
            extract_barcodes(bam_file, left_region, left_barcodes, barcodes_occ)
            left_barcodes.seek(0)
            left = set(left_barcodes.read().splitlines())
    
        #Obtain the right barcodes and store the elements in a set
        right_region = right_scaffold.chunk(args.chunk)
        right_barcodes_file = "{}{}.c{}.right.barcodes".format(right_scaffold.name, right_scaffold.orient, args.chunk)

        print("Barcodes from right chunk ({}): {}".format(right_region, right_barcodes_file))
        with open(right_barcodes_file, "w+") as right_barcodes:
            extract_barcodes(bam_file, right_region, right_barcodes, barcodes_occ)
            right_barcodes.seek(0)
            right = set(right_barcodes.read().splitlines())

        #Calculate the union 
        union_barcodes_file = "{}.{}.g{}.c{}.bxu".format(gfa_name, str(gap_label), gap.length, args.chunk)
        print("Barcodes from the union (all barcodes): " + union_barcodes_file)
        with open(union_barcodes_file, "w") as union_barcodes:
            union = left | right
            #filter barcodes by freq
            for (barcode, occurences) in barcodes_occ.items():
                if occurences >= args.freq:
                    union_barcodes.write(barcode + "\n")
        
        #----------------------------------------------------
        # GetReads
        #----------------------------------------------------
        #Union: extract the reads associated with the barcodes
        union_reads_file = "{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, str(gap_label), gap.length, args.chunk)
        print("Extracting reads associated with the barcodes from union: " + union_reads_file)
        with open(union_reads_file, "w") as union_reads:
            get_reads(reads_file, index_file, union_barcodes_file, union_reads)

        #----------------------------------------------------
        # Summary of union (barcodes and reads)
        #----------------------------------------------------
        bxu = sum(1 for line in open(union_barcodes_file, "r"))
        rbxu = sum(1 for line in open(union_reads_file, "r"))/4
        union_summary = [gap.id, gap.left, gap.right, gap.length, args.chunk, bxu, rbxu]
        
        if os.path.exists(outDir + "/{}.union.sum".format(gfa_name)):
            with open("{}.union.sum".format(gfa_name), "a") as union_sum:
                union_sum.write("\n" + '\t\t'.join(str(i) for i in union_summary))
        else:
            with open("{}.union.sum".format(gfa_name), "a") as union_sum:
                legend = ["Gap ID", "Left scaffold", "Right scaffold", "Gap size", "Chunk size", "Nb barcodes", "Nb reads"]
                union_sum.write('\t'.join(j for j in legend))
                union_sum.write("\n" + '\t\t'.join(str(i) for i in union_summary))

        #Remove the barcodes files
        subprocess.run(["rm", left_barcodes_file])
        subprocess.run(["rm", right_barcodes_file])
        subprocess.run(["rm", union_barcodes_file])

        #----------------------------------------------------
        # MindTheGap pipeline
        #----------------------------------------------------
        #Directory for saving the results from MindTheGap
        if not os.path.exists(outDir + "/mtg_results"):
            os.mkdir(outDir + "/mtg_results")
        try:
            os.chdir(outDir + "/mtg_results")
        except:
            print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
            print("Restoring the path")
            os.chdir(outDir)
        mtgDir = os.getcwd()
            
        #Execute MindTheGap fill module on the union, in breakpoint mode
        solution = False
        for k in args.k_mtg:

            #MindTheGap output directory
            os.chdir(mtgDir)
        
            #----------------------------------------------------
            # Breakpoint file, with offset of size k removed
            #----------------------------------------------------
            #variable 'ext' is the size of the extension of the gap, on both sides [by default k]
            if args.ext is None:
                ext = k
            else:
                ext = args.ext

            bkpt_file = "{}.{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, str(gap_label), gap.length, args.chunk, k)
            with open(bkpt_file, "w") as bkpt:
                line1 = ">bkpt1_GapID.{}_Gaplen.{} left_kmer.{}{}_len.{} offset_rm\n".format(str(gap_label), gap.length, left_scaffold.name, left_scaffold.orient, k)
                line2 = str(left_scaffold.sequence()[(left_scaffold.len - ext - k):(left_scaffold.len - ext)])
                line3 = "\n>bkpt1_GapID.{}_Gaplen.{} right_kmer.{}{}_len.{} offset_rm\n".format(str(gap_label), gap.length, right_scaffold.name, right_scaffold.orient, k)
                line4 = str(right_scaffold.sequence()[ext:(ext + k)])
                line5 = "\n>bkpt2_GapID.{}_Gaplen.{} left_kmer.{}{}_len.{} offset_rm\n".format(str(gap_label), gap.length, right_scaffold.name, gfapy.invert(right_scaffold.orient), k)
                line6 = str(rc(right_scaffold.sequence())[(right_scaffold.len - ext - k):(right_scaffold.len - ext)])
                line7 = "\n>bkpt2_GapID.{}_Gaplen.{} right_kmer.{}{}_len.{} offset_rm\n".format(str(gap_label), gap.length, left_scaffold.name, gfapy.invert(left_scaffold.orient), k)
                line8 = str(rc(left_scaffold.sequence())[ext:(ext + k)])
                bkpt.writelines([line1, line2, line3, line4, line5, line6, line7, line8])
            print("\nBreakpoint file (with offset of size k removed): " + bkpt_file)

            #----------------------------------------------------
            # Gapfilling
            #----------------------------------------------------
            for a in args.a_mtg:

                print("\nGapfilling of {}.{}.g{}.c{} for k={} and a={} (union)".format(gfa_name, str(gap_label), gap.length, args.chunk, k, a))
                input_file = os.path.join(outDir, union_reads_file)
                output = "{}.{}.g{}.c{}.k{}.a{}.bxu".format(gfa_name, str(gap_label), gap.length, args.chunk, k, a)
                max_nodes = args.max_nodes_mtg
                max_length = args.max_length_mtg
                if max_length == 10000 and gap.length >= 10000:
                    max_length = gap.length + 1000
                nb_cores = args.nb_cores_mtg
                max_memory = args.max_memory_mtg
                verbose = args.verbose_mtg
                mtg_fill(input_file, bkpt_file, k, a, max_nodes, max_length, nb_cores, max_memory, verbose, output)

                if os.path.getsize(mtgDir +"/"+ output + ".insertions.fasta") > 0:
                    insertion_file = os.path.abspath(mtgDir +"/"+ output + ".insertions.fasta")

                    #Modify the 'insertion_file' and save it to a new file ('input_file') so that the 'solution x/y' part appears in record.id (and not just in record.description)
                    input_file = os.path.abspath(mtgDir +"/"+ output + "..insertions.fasta")
                    with open(insertion_file, "r") as original, open(input_file, "w") as corrected:
                        records = SeqIO.parse(original, "fasta")
                        for record in records:
                            if "solution" in record.description:
                                record.id = record.id + "_sol_" + record.description.split (" ")[-1]
                            else:
                                record.id = record.id + "_sol_1/1"
                            SeqIO.write(record, corrected, "fasta")

                    #----------------------------------------------------
                    # Stats of the alignments query_seq vs reference_seq
                    #----------------------------------------------------
                    #Get the reference sequence file
                    if args.refDir is not None or args.scaffs is not None:
                        print("\nStatistical analysis...")

                        if args.refDir is not None:
                            ref_file = refDir +"/"+ str(gap_label) +".g"+ str(gap.length) + ".ingap.fasta"
                        else:
                            ref_file = scaffs_file

                        if not os.path.isfile(ref_file):
                            print("Something wrong with the specified reference file. Exception-", sys.exc_info())

                        #Do statistics on the alignments of query_seq (found gapfill seq) vs reference_seq
                        else:
                            prefix = "{}.k{}.a{}".format(str(gap_label), k, a) 
                            stats_align(input_file, ref_file, str(ext), prefix, statsDir)

                        #----------------------------------------------------
                        # Estimate quality of gapfilled sequence
                        #----------------------------------------------------
                        #Reader for alignment stats' files
                        ref_qry_file = statsDir + "/" + prefix + ".ref_qry.alignment.stats"
                        qry_qry_file = statsDir + "/" + prefix + ".qry_qry.alignment.stats"

                        if not os.path.exists(ref_qry_file) or not os.path.exists(qry_qry_file):
                            solution = False

                        else:
                            ref_qry_output = open(ref_qry_file)
                            qry_qry_output = open(qry_qry_file)

                            reader_ext_stats = csv.DictReader(ref_qry_output, \
                                                            fieldnames=("Gap", "Len_gap", "Chunk", "k", "a", "Strand", "Solution", "Len_Q", "Ref", "Len_R", \
                                                                        "Start_ref", "End_ref", "Start_qry", "End_qry", "Len_alignR", "Len_alignQ", "%_Id", "%_CovR", "%_CovQ", "Frame_R", "Frame_Q", "Quality"), \
                                                            delimiter='\t')

                            reader_revcomp_stats = csv.DictReader(qry_qry_output, \
                                                                fieldnames=("Gap", "Len_gap", "Chunk", "k", "a", "Solution1", "Len_Q1", "Solution2", "Len_Q2", \
                                                                            "Start_Q1", "End_Q1", "Start_Q2", "End_Q2", "Len_align_Q1", "Len_align_Q2", "%_Id", "%_Cov_Q1", "%_Cov_Q2", "Frame_Q1", "Frame_Q2", "Quality"), \
                                                                delimiter='\t')
                            
                            #Obtain a quality score for each gapfilled seq
                            insertion_quality_file = os.path.abspath(mtgDir +"/"+ output + ".insertions_quality.fasta")
                            with open(input_file, "r") as query, open(insertion_quality_file, "w") as qualified:
                                for record in SeqIO.parse(query, "fasta"):

                                    #quality score for stats about the extension
                                    quality_ext_left = []
                                    quality_ext_right = []
                                    for row in reader_ext_stats:
                                        if (row["Solution"] in record.id) and (("bkpt1" in record.id and row["Strand"] == "fwd") or ("bkpt2" in record.id and row["Strand"] == "rev")) and (row["Ref"] == left_scaffold.name):
                                            quality_ext_left.append(row["Quality"])
                                        elif (row["Solution"] in record.id) and (("bkpt1" in record.id and row["Strand"] == "fwd") or ("bkpt2" in record.id and row["Strand"] == "rev")) and (row["Ref"] == right_scaffold.name):
                                            quality_ext_right.append(row["Quality"])
                                    if quality_ext_left == []:
                                        quality_ext_left.append('D')
                                    if quality_ext_right == []:
                                        quality_ext_right.append('D')

                                    ref_qry_output.seek(0)

                                    #quality score for stats about the reverse complement strand
                                    quality_revcomp = []
                                    for row in reader_revcomp_stats:
                                        if ((record.id.split('_')[-1] in row["Solution1"]) and (("bkpt1" in record.id and "fwd" in row["Solution1"]) or ("bkpt2" in record.id and "rev" in row["Solution1"]))) \
                                            or ((record.id.split('_')[-1] in row["Solution2"]) and (("bkpt1" in record.id and "fwd" in row["Solution2"]) or ("bkpt2" in record.id and "rev" in row["Solution2"]))):
                                            quality_revcomp.append(row["Quality"])
                                    if quality_revcomp == []:
                                        quality_revcomp.append('D')
                                    qry_qry_output.seek(0)

                                    #global quality score
                                    quality_gapfilled_seq = min(quality_ext_left) + min(quality_ext_right) + min(quality_revcomp)
                                    
                                    record.description = "Quality " + str(quality_gapfilled_seq)
                                    SeqIO.write(record, qualified, "fasta")

                                    #If at least one good solution amongst all solution found, stop searching
                                    if re.match('^.*Quality A[AB]{2}$', record.description) or re.match('^.*Quality BA[AB]$', record.description):
                                        solution = True

                                qualified.seek(0)

                            #remove the 'input_file' once done with it
                            subprocess.run(["rm", input_file])

                            #remplace the 'insertion_file' by the 'insertion_quality_file' (which is then renamed 'insertion_file')
                            subprocess.run(["rm", insertion_file])
                            subprocess.run(['mv', insertion_quality_file, insertion_file])


                    #-------------------------------------------------------------------
                    # GFA output: case gap, solution found (=query), length(gap) known
                    #-------------------------------------------------------------------
                    #If length(gap) provided, check that length(solution found) = length(gap) +/- 10% 
                    if gap.length != 0 and solution == True:
                        with open(insertion_file, "r") as query:
                            for record in SeqIO.parse(query, "fasta"): #x records loops (x = nb of query (e.g. nb of inserted seq))
                                seq = record.seq
                                len_seq = len(seq) - 2*ext

                                if len_seq >= 0.85*gap.length and len_seq <= 1.15*gap.length:
                                    solution = True

                                    #Update GFA with only the good solutions (the ones having a good quality score)
                                    if re.match('^.*Quality A[AB]{2}$', record.description) or re.match('^.*Quality BA[AB]$', record.description):
                                        os.chdir(outDir)
                                        print("\nCreating or appending the output GFA file...")
                                        output_gfa_with_solution(outDir, record, ext, k, gap.left, gap.right, left_scaffold, right_scaffold, gfa_name, out_gfa_file)
                                    
                                    break

                                else:
                                    solution = False
    
                    #-------------------------------------------------------------------
                    # GFA output: case gap, solution found (=query), length(gap) unknown
                    #-------------------------------------------------------------------
                    #Update GFA with only the good solutions (the ones having a good quality score)
                    elif solution == True:
                        with open(insertion_file, "r") as query:
                            for record in SeqIO.parse(query, "fasta"):  #x records loops (x = nb of query (e.g. nb of inserted seq))
                                if re.match('^.*Quality A[AB]{2}$', record.description) or re.match('^.*Quality BA[AB]$', record.description):
                                    os.chdir(outDir)
                                    print("\nCreating or appending the output GFA file...")
                                    output_gfa_with_solution(outDir, record, ext, k, gap.left, gap.right, left_scaffold, right_scaffold, gfa_name, out_gfa_file)
                        
                        break


                #If no solution found, remove the 'xxx.insertions.fasta' and 'xxx.insertions.vcf' file
                else:
                    insertion_fasta = os.path.abspath(mtgDir +"/"+ output + ".insertions.fasta")
                    insertion_vcf = os.path.abspath(mtgDir +"/"+ output + ".insertions.vcf")
                    subprocess.run(["rm", insertion_fasta])
                    subprocess.run(["rm", insertion_vcf])

            if solution == True and not args.force:
                break

            #----------------------------------------------------
            # GFA output: case gap, no solution
            #----------------------------------------------------
            elif k == min(args.k_mtg) and a == min(args.a_mtg):
                #GFA output directory
                os.chdir(outDir)
                print("\nCreating or appending the output GFA file...")

                #Rewrite the current G line into GFA output
                with open("tmp.gap", "r") as tmp_gap, open(out_gfa_file, "a") as f:
                    out_gfa = gfapy.Gfa.from_file(out_gfa_file)
                    for line in tmp_gap.readlines():
                        out_gfa.add_line(line)
                    out_gfa.to_file(out_gfa_file)
        

        #Remove the tmp.gap file
        tmp_gap_file = os.path.join(outDir, "tmp.gap")
        subprocess.run(["rm", tmp_gap_file])


    #Give the output GFA file and the file containing the gapfill seq
    print("GFA file: " + out_gfa_file)


except Exception as e:
    print("\nException-")
    print(e)
    sys.exit(1)


print("\nSummary of the union: " +gfa_name+".union.sum")
print("The results from MindTheGap are saved in " + mtgDir)
print("The statistics from MTG10X are saved in " + statsDir)
print("The GFA output file and the sequences file are saved in " + outDir)