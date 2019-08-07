#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import argparse
import re
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO
import subprocess


def main():

    #----------------------------------------------------
    # Arg parser
    #----------------------------------------------------
    parser = argparse.ArgumentParser(prog="mtg10x.py", usage="%(prog)s -gfa <GFA_file> -c <chunk_size> -bam <BAM_file> -reads <reads_file> -index <index_file> [options]", \
                                     formatter_class=argparse.RawTextHelpFormatter, \
                                     description=(''' \
                                     Gapfilling with 10X data, using MindTheGap in 'breakpoint' mode
                                     ----------------------------------------------------------------
                                     To be able to execute this script, you need to install a virtual environment containing Samtools, Biopython and Gfapy.
                                     You also need to install Con10X (especially BamExtractor and GetReads) and MindTheGap 
                                     Use MTG in breakpoint mode, but taking an offset of size k
                                    
                                     [Main options]: '-gfa' : input GFA file containing the contigs paths
                                                     '-c': chunk size 
                                                     '-bam' : BAM file
                                                     '-reads': file of indexed reads
                                                     '-index': barcodes index file
                                                     '-out' : output directory [default './mtg10x_results']

                                     [MindTheGap options]: '-bkpt': breakpoint file (with possibly offset of size k removed)              
                                                           '-kmer-size' : size of a kmer [default '[51, 41, 31, 21]']
                                                           '-abundance-min' : minimal abundance threshold for solid kmers [default '[3, 2]']
                                                           '-max-nodes' : maximum number of nodes in contig graph [default '1000']
                                                           '-max-length': maximum length of gapfilling (nt) [default '10000']
                                                           '-nb-cores' : number of cores [default '4']
                                                           '-max-memory' : max memory for graph building (in MBytes) [default '8000']
                                                           '-verbose': verbosity level [default '0']
                                     '''))

    parserMain = parser.add_argument_group("[Main options]")
    parserMtg = parser.add_argument_group("[MindTheGap option]")

    parserMain.add_argument('-gfa', action="store", dest="gfa", help="input GFA file containing the contigs paths (format: xxx.gfa)", required=True)
    parserMain.add_argument('-c', action="store", dest="chunk", type=int, help="size of the chunk for gapfilling", required=True)
    parserMain.add_argument('-bam', action="store", dest="bam", help="BAM file containing the reads of the individual", required=True)
    parserMain.add_argument('-reads', action="store", dest="reads", help="file of indexed reads", required=True)
    parserMain.add_argument('-index', action="store", dest="index", help="barcodes index file", required=True)
    parserMain.add_argument('-out', action="store", dest="out_dir", default="./mtg10x_results", help="output directory for result files")

    parserMtg.add_argument('-bkpt', action="store", dest="bkpt", help="breakpoint file in fasta format")
    parserMtg.add_argument('-kmer-size', action="store", dest= "k_mtg", default=[51, 41, 31, 21],  nargs='*', type=int, help="kmer size used for gapfilling")
    parserMtg.add_argument('-abundance-min', action="store", dest="a_mtg", default=[3, 2], nargs='*', type=int, help="minimal abundance of kmers used for gapfilling")
    parserMtg.add_argument('-max-nodes', action="store", dest="max_nodes_mtg", type=int, default=1000, help="maximum number of nodes in contig graph")
    parserMtg.add_argument('-max-length', action="store", dest="max_length_mtg", type=int, default=10000, help="maximum length of gapfilling (nt)")
    parserMtg.add_argument('-nb-cores', action="store", dest="nb_cores_mtg", type=int, default=4, help="number of cores")
    parserMtg.add_argument('-max-memory', action="store", dest="max_memory_mtg", type=int, default=8000, help="max memory for graph building (in MBytes)")

    global args
    args = parser.parse_args()

    if re.match('^.*\.gfa$', args.gfa) is None:
        parser.error("The suffix of the GFA file should be: '.gfa'")
    elif not os.path.exists(args.gfa):
        parser.error("The path of the GFA file doesn't exist")

    if re.match('^.*\.bam$', args.bam) is None:
        parser.error("The suffix of the BAM file should be: '.bam'")
    elif not os.path.exists(args.bam): 
        parser.error("The path of the BAM file doesn't exist")

    if not os.path.exists(args.reads):
        parser.error("The path of the file of indexed reads doesn't exist")
    
    if not os.path.exists(args.index):
        parser.error("The path of the barcodes index file doesn't exist")

    #----------------------------------------------------
    # Input files
    #----------------------------------------------------
    gfa_file = args.gfa
    gfa_name = gfa_file.split('/')[-1]
    bam_file = args.bam
    print("\nInput GFA file: " + gfa_file)
    print("BAM file: " + bam_file)
    print("File of indexed reads: " + args.reads)
    print("Barcodes index file: " + args.index)

    #----------------------------------------------------
    # Directory for saving results
    #----------------------------------------------------
    cwd = os.getcwd() 
    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    try:
        os.chdir(args.out_dir)
    except:
        print("Something wrong with specified directory. Exception-", sys.exc_info())
        print("Restoring the path")
        os.chdir(cwd)
    cwd = os.getcwd()
    print("\nThe results are saved in " + cwd)

    #----------------------------------------------------
    # Gapfilling pipeline
    #----------------------------------------------------
    try:
        gfa = gfapy.Gfa.from_file(gfa_file)
        for gap in gfa.gaps:
            global gap_len
            l_contig = gap.sid1
            r_contig = gap.sid2
            gap_len = gap.disp
            gap_id = gap.gid
            
            if gap_id != '*':
                print("\nWORKING ON: {}.{}.g{}.c{}".format(gfa_name, gap_id, gap_len, args.chunk))
            elif gap_id == '*':
                print("\nWORKING ON: {}.{}-{}.g{}.c{}".format(gfa_name, l_contig, r_contig, gap_len, args.chunk))
            
            #----------------------------------------------------
            # GFA handling
            #----------------------------------------------------
            #Obtain the name, the sequence and the sequence length of the left contig
            (l_name, l_len, l_seq, l_orient) = gfa_handle(l_contig)

            #Obtain the name, the sequence and the sequence length of the right contig
            (r_name, r_len, r_seq, r_orient) = gfa_handle(r_contig)

            #Obtain the left region and right region on which to extract the barcodes
            if args.chunk > l_len or args.chunk > r_len:
                sys.stderr.write("\nError: The chunk size must be smaller than the sequence length of the contigs")
                sys.exit(2)
            else:
                l_start = l_len - args.chunk 
                l_end = l_len
                r_start = 0 
                r_end = args.chunk
           
            #----------------------------------------------------
            # BamExtractor
            #----------------------------------------------------
            #Obtain the left barcodes and store the elements in a set
            left_region = "{}:{}-{}".format(l_name, l_start, l_end)
            print("\nBarcodes from left chunk ({}): {}_{}.g{}.c{}.left.barcodes".format(left_region, l_name, l_orient, gap_len, args.chunk))
            with open("{}_{}.g{}.c{}.left.barcodes".format(l_name, l_orient, gap_len, args.chunk), "w+") as left_barcodes:
                bam_extract(bam_file, left_region, left_barcodes)
                left_barcodes.seek(0)
                l_barcodes = left_barcodes.read()
                l = set(l_barcodes.splitlines())

            #Obtain the right barcodes and store the elements in a set
            right_region = "{}:{}-{}".format(r_name, r_start, r_end)
            print("Barcodes from right chunk ({}): {}_{}.g{}.c{}.right.barcodes".format(right_region, r_name, r_orient, gap_len, args.chunk))
            with open("{}_{}.g{}.c{}.right.barcodes".format(r_name, r_orient, gap_len, args.chunk), "w+") as right_barcodes:
                bam_extract(bam_file, right_region, right_barcodes)
                right_barcodes.seek(0)
                r_barcodes = right_barcodes.read()
                r = set(r_barcodes.splitlines())

            #Calculate the union
            if gap_id != '*':
                print("Barcodes from the union (all barcodes): {}.{}.g{}.c{}.bxu".format(gfa_name, gap_id, gap_len, args.chunk))
                with open("{}.{}.g{}.c{}.bxu".format(gfa_name, gap_id, gap_len, args.chunk), "w") as union_barcodes:    
                    union = l | r
                    union_barcodes.write('\n'.join(seq for seq in union))

            elif gap_id == '*':
                print("Barcodes from the union (all barcodes): {}.{}-{}.g{}.c{}.bxu".format(gfa_name, l_contig, r_contig, gap_len, args.chunk))
                with open("{}.{}-{}.g{}.c{}.bxu".format(gfa_name, l_contig, r_contig, gap_len, args.chunk), "w") as union_barcodes:    
                    union = l | r
                    union_barcodes.write('\n'.join(seq for seq in union))

            #----------------------------------------------------
            # GetReads
            #----------------------------------------------------
            #Union: extract the reads associated with the barcodes
            if gap_id != '*':
                print("Extracting reads associated with the barcodes from union: {}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, gap_id, gap_len, args.chunk))
                u_barcodes = cwd + "/{}.{}.g{}.c{}.bxu".format(gfa_name, gap_id, gap_len, args.chunk)
                with open("{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, gap_id, gap_len, args.chunk), "w") as union_reads:
                    get_reads(u_barcodes, union_reads)

            elif gap_id == '*':
                print("Extracting reads associated with the barcodes from union: {}.{}-{}.g{}.c{}.rbxu.fastq".format(gfa_name, l_contig, r_contig, gap_len, args.chunk))
                u_barcodes = cwd + "/{}.{}-{}.g{}.c{}.bxu".format(gfa_name, l_contig, r_contig, gap_len, args.chunk)
                with open("{}.{}-{}.g{}.c{}.rbxu.fastq".format(gfa_name, l_contig, r_contig, gap_len, args.chunk), "w") as union_reads:
                    get_reads(u_barcodes, union_reads)

            #----------------------------------------------------
            # Summary of union (barcodes and reads)
            #----------------------------------------------------
            if gap_id != '*':
                bxu = sum(1 for line in open("{}.{}.g{}.c{}.bxu".format(gfa_name, gap_id, gap_len, args.chunk), "r"))
                rbxu = sum(1 for line in open("{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, gap_id, gap_len, args.chunk), "r"))/4

            elif gap_id == '*':
                bxu = sum(1 for line in open("{}.{}-{}.g{}.c{}.bxu".format(gfa_name, l_contig, r_contig, gap_len, args.chunk), "r"))
                rbxu = sum(1 for line in open("{}.{}-{}.g{}.c{}.rbxu.fastq".format(gfa_name, l_contig, r_contig, gap_len, args.chunk), "r"))/4

            u_summary = [gap_id, l_contig, r_contig, gap_len, args.chunk, bxu, rbxu]
            if os.path.exists(cwd + "/{}.union.sum".format(gfa_name)):
                with open("{}.union.sum".format(gfa_name), "a") as union_sum:
                    union_sum.write("\n" + '\t\t'.join(str(i) for i in u_summary))
            else:
                with open("{}.union.sum".format(gfa_name), "a") as union_sum:
                    legend = ["Gap ID", "Left contig", "Right contig", "Gap size", "Chunk size", "Nb barcodes", "Nb reads associated"]
                    union_sum.write('\t'.join(j for j in legend))
                    union_sum.write("\n" + '\t\t'.join(str(i) for i in u_summary))

            #----------------------------------------------------
            # MindTheGap pipeline
            #----------------------------------------------------
            #Directory for saving the results from MindTheGap
            if not os.path.exists(cwd + "/mtg_results"):
                os.mkdir(cwd + "/mtg_results")
            try:
                os.chdir(cwd + "/mtg_results")
            except:
                print("\nSomething wrong with specified directory. Exception-", sys.exc_info())
                print("Restoring the path")
                os.chdir(cwd)
            pwd = os.getcwd()
             
            #Execute MindTheGap fill module on the union, in breakpoint mode
            solution = False
            for k in args.k_mtg:
            
                #----------------------------------------------------
                # Breakpoint file, with offset of size k removed
                #----------------------------------------------------
                if gap_id != '*':
                    with open("{}.{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, gap_id, gap_len, args.chunk, k), "w") as bkpt:
                        line1 = ">bkpt1_GapID.{}_Gaplen.{} left_kmer.{}{}_len.{} offset_rm\n".format(gap_id, gap_len, l_contig, l_contig.orient, k)
                        line2 = str(l_seq[(l_len - 2*k):(l_len - k)])
                        line3 = "\n>bkpt1_GapID.{}_Gaplen.{} right_kmer.{}{}_len.{} offset_rm\n".format(gap_id, gap_len, r_contig, r_contig.orient, k)
       	       	        line4 = str(r_seq[k:2*k])
                        line5 = "\n>bkpt2_GapID.{}_Gaplen.{} left_kmer.{}{}_len.{} offset_rm\n".format(gap_id, gap_len, r_contig, gfapy.invert(r_contig.orient), k)
                        line6 = str(rc(r_seq)[(r_len - 2*k):(r_len - k)])
                        line7 = "\n>bkpt2_GapID.{}_Gaplen.{} right_kmer.{}{}_len.{} offset_rm\n".format(gap_id, gap_len, l_contig, gfapy.invert(l_contig.orient), k)
                        line8 = str(rc(l_seq)[k:2*k])
       	       	        bkpt.writelines([line1, line2, line3, line4, line5, line6, line7, line8])
                    bkpt_file = pwd + "/{}.{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, gap_id, gap_len, args.chunk, k)
                
                elif gap_id == '*':
                    with open("{}.{}-{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, l_contig, r_contig, gap_len, args.chunk, k), "w") as bkpt:
                        line1 = ">bkpt1_Gaplen.{0}_Contigs.{1}-{2} left_kmer.{1}{3}_len.{4} offset_rm\n".format(gap_len, l_contig, r_contig, l_contig.orient, k)
                        line2 = str(l_seq[(l_len - 2*k):(l_len - k)])
                        line3 = "\n>bkpt1_Gaplen.{0}_Contigs.{1}-{2} right_kmer.{2}{3}_len.{4} offset_rm\n".format(gap_len, l_contig, r_contig, r_contig.orient, k)
       	       	        line4 = str(r_seq[k:2*k])
                        line5 = "\n>bkpt2_Gaplen.{0}_Contigs.{1}-{2} left_kmer.{2}{3}_len.{4} offset_rm\n".format(gap_len, l_contig, r_contig, gfapy.invert(r_contig.orient), k)
                        line6 = str(rc(r_seq)[(r_len - 2*k):(r_len - k)])
                        line7 = "\n>bkpt2_Gaplen.{0}_Contigs.{1}-{2} right_kmer.{1}{3}_len.{4} offset_rm\n".format(gap_len, l_contig, r_contig, gfapy.invert(l_contig.orient), k)
                        line8 = str(rc(l_seq)[k:2*k])
       	       	        bkpt.writelines([line1, line2, line3, line4, line5, line6, line7, line8])
                    bkpt_file = pwd + "/{}.{}-{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, l_contig, r_contig, gap_len, args.chunk, k)

                print("\nBreakpoint file (with offset of size k removed): " + bkpt_file)

                #----------------------------------------------------
                # Gapfilling
                #----------------------------------------------------
                for a in args.a_mtg:

                    if gap_id != '*':
                        print("\nGapfilling of {}.{}.g{}.c{} for k={} and a={} (union)".format(gfa_name, gap_id, gap_len, args.chunk, k, a))
                        input_file = cwd + "/{}.{}.g{}.c{}.rbxu.fastq".format(gfa_name, gap_id, gap_len, args.chunk)
                        bkpt_file = pwd + "/{}.{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, gap_id, gap_len, args.chunk, k)
                        output = "{}.{}.g{}.c{}.k{}.a{}.bxu".format(gfa_name, gap_id, gap_len, args.chunk, k, a)
                        mtg_gapfill(input_file, bkpt_file, k, a, output)

                        if os.path.getsize(pwd + "/{}.{}.g{}.c{}.k{}.a{}.bxu.insertions.fasta".format(gfa_name, gap_id, gap_len, args.chunk, k, a)) > 0:
                            solution = True
                            break

                    elif gap_id == '*':
                        print("\nGapfilling of {}.{}-{}.g{}.c{} for k={} and a={} (union)".format(gfa_name, l_contig, r_contig, gap_len, args.chunk, k, a))
                        input_file = cwd + "/{}.{}-{}.g{}.c{}.rbxu.fastq".format(gfa_name, l_contig, r_contig, gap_len, args.chunk)
                        bkpt_file = pwd + "/{}.{}-{}.g{}.c{}.k{}.offset_rm.bkpt.fasta".format(gfa_name, l_contig, r_contig, gap_len, args.chunk, k)
                        output = "{}.{}-{}.g{}.c{}.k{}.a{}.bxu".format(gfa_name, l_contig, r_contig, gap_len, args.chunk, k, a)
                        mtg_gapfill(input_file, bkpt_file, k, a, output)

                        if os.path.getsize(pwd + "/{}.{}-{}.g{}.c{}.k{}.a{}.bxu.insertions.fasta".format(gfa_name, l_contig, r_contig, gap_len, args.chunk, k, a)) > 0:
                            solution = True
                            break
                
                if solution == True:
                    break
            
            if re.match('^.*\/mtg_results$', pwd) is None:
                os.chdir(cwd)
            else:
                os.chdir(os.path.join(pwd, '..'))


    except Exception as e:
        print("\nException-")
        print(e)
        sys.exit(1)

    
    print("\nSummary of the union: " +gfa_name+".union.sum")
    print("The results from MindTheGap are saved in " + pwd)


#----------------------------------------------------
# gfa_handle function
#----------------------------------------------------              
def gfa_handle(contig):
    name_ = contig.name
    len_ = contig.line.slen
    seq_path = contig.line.UR

    for record in SeqIO.parse(seq_path, "fasta"): 
        if contig.orient == "+":
            sequence = record.seq
            orient = "fwd"
        elif contig.orient == "-":
            sequence = rc(record.seq)
            orient = "rev"
    return name_, len_, sequence, orient

#----------------------------------------------------
# bam_extract function
#----------------------------------------------------
#Function to bam_extract the barcodes from the chunks with BamExtractor 
#(remove the '-1' at the end of the sequences (p2), and keep only the barcodes observed more than once (p3, p4, p5))
def bam_extract(bam, region, barcodes):
    command = ["BamExtractor", bam, region]
    p1 = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2 = subprocess.Popen(["cut", "-d", "-", "-f1"], stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(["sort"], stdin=p2.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p2.stdout.close()
    p4 = subprocess.Popen(["uniq", "-c"], stdin=p3.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p3.stdout.close()
    p5 = subprocess.Popen(["awk", "$1 > 1 {print $2}"], stdin=p4.stdout, stdout=barcodes, stderr=barcodes)
    p4.stdout.close()
    p5.wait()
    return barcodes

#----------------------------------------------------
# get_reads function
#----------------------------------------------------
#Function to extract the reads associated with the barcodes
def get_reads(barcodes, out_reads):    
    reads = args.reads
    index = args.index
    command = ["GetReads", "-reads", reads, "-index", index, "-barcodes", barcodes]
    subprocess.run(command, stdout=out_reads, stderr=out_reads)
    return out_reads

#----------------------------------------------------
# mtg_gapfill function
#----------------------------------------------------
#Function to execute MindTheGap fill module
def mtg_gapfill(input_file, bkpt, k, a, output_prefix):
    max_nodes = args.max_nodes_mtg
    max_length = args.max_length_mtg
    if max_length == 10000 and gap_len >= 10000:
        max_length = gap_len + 1000
    nb_cores = args.nb_cores_mtg
    max_memory = args.max_memory_mtg
    command = ["MindTheGap", "fill", "-in", input_file, "-bkpt", bkpt, "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(max_nodes), "-max-length", str(max_length), \
               "-nb-cores", str(nb_cores), "-max-memory", str(max_memory), "-out", output_prefix]
    subprocess.run(command)
    output = subprocess.check_output(command)
    subprocess.run("rm -f *.h5", shell=True)
    return output
    #!!Non: enregistrer l'output de sortie dans un fichier plutôt, fichier log ??


if __name__ == "__main__":
    main()