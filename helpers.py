#!/usr/bin/env python3
#*****************************************************************************
#  Name: MTG-Link
#  Description: gap-filling tool for draft genome assemblies, dedicated to 
#  linked read data generated by 10XGenomics Chromium technology.
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

import os
import sys
import re
import subprocess
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO
from datetime import datetime


#----------------------------------------------------
# Gap class
#----------------------------------------------------
class Gap:
    '''
    Class defining a gap characterized by:
    - its ID
    - its length
    - its left flanking sequence's name
    - its right flanking sequence's name
    '''

    #Constructor
    def __init__(self, gap):
        self._identity = gap.gid
        self._length = gap.disp
        self._left = gap.sid1
        self._right = gap.sid2

    #Accessors
    def _get_identity(self):
        '''Method to be call when we want to access the attribute "identity"'''
        return self._identity
    def _get_length(self):
        '''Method to be call when we want to access the attribute "length"'''
        return self._length
    def _get_left(self):
        '''Method to be call when we want to access the attribute "left"'''
        return self._left
    def _get_right(self):
        '''Method to be call when we want to access the attribute "right"'''
        return self._right

    #Properties
    identity = property(_get_identity)
    length = property(_get_length)
    left = property(_get_left)
    right = property(_get_right)

    #Method "__getattr__"
    def __getattr__(self, attr):
        '''If Python doesn't find the attribute "attr", it calls this method and print an alert'''
        print("WARNING: There is no attribute {} here !".format(attr))

    #Method "__delattr_"
    def __delattr_(self, attr):
        '''We can't delete an attribute, we raise the exception AttributeError'''
        raise AttributeError("You can't delete attributes from this class")
    
    #Method "label"
    def label(self):
        '''Method to label the gap'''
        if self._identity == "*":
            return str(self.left) +"_"+ str(self.right)
        else:
            return str(self.identity)

    #Method "info"
    def info(self):
        '''Method to get some information on the gap'''
        if self.identity == "*":
            print("WORKING ON GAP: between contigs {} & {}; length {}\n".format(self.left, self.right, self.length))
        else:
            print("WORKING ON GAP: {}; length {}\n".format(self.identity, self.length))

    #Method "__repr__"
    def __repr__(self):
        return "Gap: id ({}), length ({}), left flanking seq ({}), right flanking seq ({})".format(self.identity, self.length, self.left, self.right)


#----------------------------------------------------
# Scaffold class
#----------------------------------------------------
class Scaffold(Gap):
    '''
    Class defining a scaffold characterized by:
    - the gap it is linked to
    - its name
    - its orientation
    - its length
    - the path of its sequence
    '''

    #Constructor
    def __init__(self, gap, scaffold, gfa_file):
        super().__init__(gap)
        self.gap = gap
        self.scaffold = scaffold
        self._name = scaffold.name
        self._orient = scaffold.orient
        self._slen = scaffold.line.slen
        self._seq_path = scaffold.line.UR
        self.gfa_file = gfa_file
    
    #Accessors
    def _get_name(self):
        '''Method to be call when we want to access the attribute "name"'''
        return self._name
    def _get_orient(self):
        '''Method to be call when we want to access the attribute "orient"'''
        return self._orient
    def _get_slen(self):
        '''Method to be call when we want to access the attribute "slen"'''
        return self._slen
    def _get_seq_path(self):
        '''Method to be call when we want to access the attribute "seq_path"'''
        return self._seq_path

    #Properties
    name = property(_get_name)
    orient = property(_get_orient)
    slen = property(_get_slen)
    seq_path = property(_get_seq_path)

    #Method "__getattr__"
    def __getattr__(self, attr):
        '''If Python doesn't find the attribute "attr", it calls this method and print an alert'''
        print("WARNING: There is no attribute {} here !".format(attr))

    #Method "__delattr_"
    def __delattr_(self, attr):
        '''We can't delete an attribute, we raise the exception AttributeError'''
        raise AttributeError("You can't delete attributes from this class")

    #Method "sequence"
    def sequence(self):
        '''Method to get the sequence of the scaffold'''
        #if relative path
        if not str(self.seq_path).startswith('/'):
            seq_link = str('/'.join(str(self.gfa_file).split('/')[:-1])) +"/"+ str(self.seq_path)
        #if absolute path
        else:
            seq_link = self.seq_path
            
        #get the sequence of the scaffold
        for record in SeqIO.parse(seq_link, "fasta"):
            if re.match(self.name, record.id):
                if self.orient == "+":
                    return record.seq
                elif self.orient == "-":
                    return rc(record.seq)

    #Method "chunk"
    def chunk(self, c):
        '''Method to get the region of the chunk'''
        #----------------------------------------------------
        # For simulated datasets
        #----------------------------------------------------
        if ('-L' in self.name) or ('-R' in self.name):
            #if left scaffold
            if self.scaffold == self.left:
                start = self._slen - c
                end = self._slen
            #if right scaffold
            elif self.scaffold == self.right:
                start = self.slen + self.length
                end = self.slen + self.length + c
            contig_name = str(self.name).split("-")[0]
            return str(contig_name) +":"+ str(start) +"-"+ str(end)
        #----------------------------------------------------
        # For real datasets
        #----------------------------------------------------
        else:
            #if left_fwd or right_rev
            if (self.orient == "+" and self.scaffold == self.left) or (self.orient == "-" and self.scaffold == self.right):
                start = self.slen - c
                end = self.slen
            #if right_fwd or left_rev
            elif (self.orient == "+" and self.scaffold == self.right) or (self.orient == "-" and self.scaffold == self.left):
                start = 0
                end = c
            return str(self.name) +":"+ str(start) +"-"+ str(end)

    #Method "__repr__"
    def __repr__(self):
        return "Scaffold: name ({}), orientation ({}), length ({}), sequence's file ({})".format(self.name, self.orient, self.slen, self.seq_path)

    


#----------------------------------------------------
# extract_barcodes function
#----------------------------------------------------
def extract_barcodes(bam, gap_label, region, barcodes_occ):
    """
    To extract the barcodes of reads mapping on chunks, with `LRez extract`:
    - it takes as input the BAM file, the gap label, the chunk region on which to extract the barcodes, and the dictionary 'barcodes_occ'
    - it outputs the updated dictionary 'barcodes_occ' containing the occurences for each barcode extracted on the chunk region
    """

    #LRez extract
    ## parameter -d to include duplicate barcodes
    command = ["LRez", "extract", "--bam", bam, "--region", region, "-d"]
    extractLog = str(gap_label) + "_extract.log"
    tmp_barcodes_file = str(gap_label) + "_extract_stdout.txt"

    with open(tmp_barcodes_file, "w+") as f, open(extractLog, "a") as log:
        subprocess.run(command, stdout=f, stderr=log)
        f.seek(0)

        #Save the barcodes and their occurences in the dict 'barcodes_occ'
        for line in f.readlines():
            #remove the '-1' at the end of the sequence
            barcode_seq = line.split('-')[0]
            #count occurences of each barcode and save them in the dict 'barcodes_occ'
            if barcode_seq in barcodes_occ:
                barcodes_occ[barcode_seq] += 1
            else:
                barcodes_occ[barcode_seq] = 1

    #remove the raw files obtained from LRez extract
    subprocess.run(["rm", tmp_barcodes_file])
    if os.path.getsize(extractLog) == 0:
        subprocess.run(["rm", extractLog])

    return barcodes_occ


#----------------------------------------------------
# get_reads function
#----------------------------------------------------
def get_reads(reads, index, gap_label, barcodes, out_reads):
    """
    To extract the the reads associated to the barcodes:
    - it takes as input the reads file, the barcodes index file, the gap label, the file containing the barcodes of the union, and the name of the output file containing the reads of the union
    - it outputs the file containing the reads of the union
    """

    #LRez query fastq
    command = ["LRez", "query", "fastq", "--fastq", reads, "--index", index, "--list", barcodes]
    getreadsLog = str(gap_label) + "_getreads.log"

    with open(getreadsLog, "a") as log:
        subprocess.run(command, stdout=out_reads, stderr=log)

    #remove the raw files obtained from LRez query fastq
    if os.path.getsize(getreadsLog) == 0:
        subprocess.run(["rm", getreadsLog])

    return out_reads


#----------------------------------------------------
# mtg_fill function
#----------------------------------------------------
'''
To execute MindTheGap fill:
    - it takes as input the gap label, the file containing the reads of the union, the bkpt file, the kmer value, the abudance threshold for solid kmers values, the maximum number of nodes in the contig graph,
      the maximum length of gap-filling (bp), the number of cores and the max memory (in MBytes) for graph building, the verbosity level and the prefix of the output file on which to save the results
    - it outputs the results obtained from MindTheGap
'''
def mtg_fill(gap_label, input_file, bkpt, k, a, max_nodes, max_length, nb_cores, max_memory, verbose, output_prefix):
    if max_memory == 0:
        command = ["/home/genouest/inra_umr1349/aguichard/Gapfilling/DBG/scripts/MindTheGap/build/bin/MindTheGap", "fill", "-in", input_file, "-bkpt", bkpt, "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(max_nodes), "-max-length", str(max_length), \
                    "-nb-cores", str(nb_cores), "-verbose", str(verbose), "-fwd-only", "-out", output_prefix]
    else:
        command = ["/home/genouest/inra_umr1349/aguichard/Gapfilling/DBG/scripts/MindTheGap/build/bin/MindTheGap", "fill", "-in", input_file, "-bkpt", bkpt, "-kmer-size", str(k), "-abundance-min", str(a), "-max-nodes", str(max_nodes), "-max-length", str(max_length), \
                    "-nb-cores", str(nb_cores), "-max-memory", str(max_memory), "-verbose", str(verbose), "-fwd-only", "-out", output_prefix]
    mtgfillLog = str(gap_label) + "_mtgfill.log"

    with open(mtgfillLog, "a") as log:
        subprocess.run(command, stderr=log)
        output = subprocess.check_output(command)

    #remove the raw files obtained from MindTheGap
    if os.path.getsize(mtgfillLog) <= 0:
        subprocess.run(["rm", mtgfillLog])

    return output


#----------------------------------------------------
# stats_align function
#----------------------------------------------------
'''
To perform statistics on the alignment between a reference sequence and query sequences
    - it takes as input the gap label, the file containing the gap-filled sequences obtained from MindTheGap, the file containing either the reference sequence or the flanking contigs' sequences, 
      the size of the extension of the gap, the prefix name of the output files, the name of the output directory for saving the results
'''
def stats_align(gap_label, qry_file, ref_file, ext, prefix, out_dir):
    scriptPath = sys.path[0]
    stats_align_command = os.path.join(scriptPath, "stats_alignment.py")
    command = [stats_align_command, "-qry", qry_file, "-ref", ref_file, "-ext", ext, "-p", prefix, "-out", out_dir]
    statsLog = str(gap_label) + "_stats_align.log"

    with open(statsLog, "a") as log:
        subprocess.run(command, stderr=log)

    #remove the raw file obtained from statistics
    if os.path.getsize(statsLog) <= 0:
        subprocess.run(["rm", statsLog])


#----------------------------------------------------
# get_position_for_edges function
#----------------------------------------------------
def get_position_for_edges(orient1, orient2, length1, length2, ext):
    #Same orientation
    if orient1 == orient2:

        #forward orientations
        if orient1 == "+":
            beg1 = str(length1 - ext)
            end1 = str(length1) + "$"  
            beg2 = str(0)
            end2 = str(ext)

        #reverse orientations
        elif orient1 == "-":
            beg1 = str(0)
            end1 = str(ext)
            beg2 = str(length2 - ext)
            end2 = str(length2) + "$"

    #Opposite orientation
    elif orient1 != orient2:

        #first seq in fwd orientation and second seq in rev orientation
        if orient1 == "+":
            beg1 = str(length1 - ext)
            end1 = str(length1) + "$"
            beg2 = str(length2 - ext)
            end2 = str(length2) + "$"

        #first seq in rev orientation and first seq in fwd orientation
        elif orient1 == "-":
            beg1 = str(0)
            end1 = str(ext)
            beg2 = str(0)
            end2 = str(ext)

    positions = [beg1, end1, beg2, end2]
    return positions


#----------------------------------------------------
# get_output_for_gfa function
#----------------------------------------------------
#Function to get the ouput variables to then output the GFA when a solution is found for a gap
def get_output_for_gfa(record, ext, k, s1, s2, left_scaffold, right_scaffold):
    seq = record.seq
    length_seq = len(seq)
    sol_name = record.id.split("_")[-1] + ".k" + str(k)
    orient_sign = "+"
    orient = "fwd"
    quality = record.description.split('Quality ')[1]

    if "bkpt2" in str(record.id):
        orient_sign = "-"
        orient = "rev"

    sol_name = str(s1) +":"+ str(s2) + "_gf" + sol_name +"_"+ orient
    solution = sol_name + orient_sign

    pos_1 = get_position_for_edges(left_scaffold.orient, orient_sign, left_scaffold.slen, length_seq, ext)
    pos_2 = get_position_for_edges(orient_sign, right_scaffold.orient, length_seq, right_scaffold.slen, ext)


    output_for_gfa = [sol_name, length_seq, str(seq), solution, pos_1, pos_2, quality]
    return output_for_gfa


#----------------------------------------------------
# update_gfa_with_solution function
#----------------------------------------------------
#Function to update the GFA when a solution is found for a gap
def update_gfa_with_solution(outDir, gfa_name, output_for_gfa, gfa_output_file):

    #Variables input
    sol_name = output_for_gfa[0]
    length_seq = output_for_gfa[1]
    seq = output_for_gfa[2]
    solution = output_for_gfa[3]
    pos_1 = output_for_gfa[4]
    pos_2 = output_for_gfa[5]
    s1 = sol_name.split(':')[0]
    s2 = (sol_name.split(':')[1]).split('_gf')[0]
    quality = output_for_gfa[6]

    print("Updating the GFA file with the solution: " + sol_name)

    #Save the found seq to a file containing all gapfill seq
    gapfill_file = gfa_name + ".gapfill_seq.fasta"
    with open(gapfill_file, "a") as seq_fasta:
        seq_fasta.write(">{} _ len_{}_qual_{} ".format(sol_name, length_seq, quality))
        seq_fasta.write("\n" + seq + "\n")

    with open(gfa_output_file, "a") as f:
        '''
        #Not used anymore, as I check in the gapfilling pipeline that the length of the solution found is > 2*ext to stop searching for a gapfilled seq
        if length_seq < 2*ext:
            print("Query length is too short (<2*ext): overlap of source and destination read")
            print("Rewriting the gap line to the GFA output file...")

            #Rewrite the current G line into GFA output
            with open("tmp.gap", "r") as tmp_gap, open(gfa_output_file, "a") as f:
                out_gfa = gfapy.Gfa.from_file(gfa_output_file)
                for line in tmp_gap.readlines():
                    out_gfa.add_line(line)
                out_gfa.to_file(gfa_output_file)

        else:
        '''
        #Add the found seq (query seq) to GFA output (S line)
        out_gfa = gfapy.Gfa.from_file(gfa_output_file)
        out_gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(sol_name, length_seq, os.path.join(outDir, gapfill_file)))

        #Write the two corresponding E lines into GFA output
        out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(s1, solution, pos_1[0], pos_1[1], pos_1[2], pos_1[3]))
        out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(solution, s2, pos_2[0], pos_2[1], pos_2[2], pos_2[3]))

        out_gfa.to_file(gfa_output_file)

        return gapfill_file