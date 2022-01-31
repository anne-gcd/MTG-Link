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

"""Module 'helpers.py': Classes and general functions

The module 'helpers.py' contains the classes and general functions used in the gap-filling pipeline MTG-Link.
"""

import os
import sys
import re
import subprocess
import gfapy
import pysam
from gfapy.sequence import rc
from Bio import SeqIO
#from datetime import datetime


#----------------------------------------------------
# Gap class
#----------------------------------------------------
class Gap:
    """
    Class defining a gap characterized by:
    - its ID
    - its length
    - its left flanking sequence's name
    - its right flanking sequence's name
    """

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
            print("WORKING ON GAP: between contigs {} & {}; length {}".format(self.left, self.right, self.length))
        else:
            print("WORKING ON GAP: {}; length {}".format(self.identity, self.length))

    #Method "__repr__"
    def __repr__(self):
        return "Gap: id ({}), length ({}), left flanking seq ({}), right flanking seq ({})".format(self.identity, self.length, self.left, self.right)


#----------------------------------------------------
# Scaffold class
#----------------------------------------------------
class Scaffold(Gap):
    """
    Class defining a scaffold characterized by:
    - the gap it is linked to
    - its name
    - its orientation
    - its length
    - the path of its sequence
    """

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
        # For gaps into scaffolds' sequences
        #----------------------------------------------------
        if ('-L' in self.name) or ('-R' in self.name):
            coordsOnScaffold = re.findall(r'[0-9]+:[0-9]+', str(self.name))[0]

            #if left scaffold
            if self.scaffold == self.left:
                start = str(coordsOnScaffold).split(':')[1] - c
                end = str(coordsOnScaffold).split(':')[1]

            #if right scaffold
            elif self.scaffold == self.right:
                start = str(coordsOnScaffold).split(':')[0]
                end = str(coordsOnScaffold).split(':')[0] + c

            #NB: The 'contig_name' should match the contig name on the BAM file
            contig_name = re.split(r'_[0-9]+:[0-9]+', str(self.name))[0]
            return str(contig_name) +":"+ str(start) +"-"+ str(end)

        #----------------------------------------------------
        # For gaps between scaffolds' sequences
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
# Graph class
#----------------------------------------------------
class Graph:
    """The class 'Graph' contains all the attributes, properties and methods to create a Graph object.

    The class 'Graph' initializes a Graph object.
    A dictionary will be used for storing the nodes and their corresponding neighbouring nodes.
    """
    # Constructor.
    def __init__(self, graph_dict):
        self._graph = graph_dict

    # Accessor.
    def _get_graph(self):
        '''Method to be call when we want to access the attribute "graph"'''
        return self._graph

    # Property.
    graph = property(_get_graph)

    # Method "__getattr__".
    def __getattr__(self, attr):
        '''If Python doesn't find the attribute "attr", it calls this method and print an alert'''
        print("WARNING: There is no attribute {} here !".format(attr))

    # Method "__delattr_".
    def __delattr_(self, attr):
        '''We can't delete an attribute, we raise the exception AttributeError'''
        raise AttributeError("You can't delete attributes from this class")

    # Method "add_node".
    def add_node(self, node):
        '''Method to add the node 'node' to the graph if it's not already in the graph'''
        if node not in self.graph:
            self.graph[node] = []

    # Method "add_edge".
    def add_edge(self, edge):
        '''Method to add an 'edge' between a node and its neighbours'''
        (source_node, neighbour, overlap) = edge
        #add an edge between the source_node and its neighbour node, with their corresponding overlap length
        if source_node in self.graph:
            self.graph[source_node].append([neighbour, overlap])
        else:
            self.graph[source_node] = [[neighbour, overlap]]

    # Method "nodes".
    def nodes(self):
        '''Method to return a list of the graph' nodes'''
        nodes = []
        for node in self.graph.keys():
            nodes.append(node)
        return nodes

    # Method "edges".
    def edges(self):
        '''Method to return a list of the graph' edges, represented as a set, with one node (a loop back to the node) or two nodes and their overlapping length'''
        edges = []
        for node in self.graph:
            for neighbour in self.graph[node]:
                if [node, neighbour[0], neighbour[1]] not in edges:
                    edges.append([node, neighbour[0], neighbour[1]])
        return edges

    # Method "create_graph_from_extensions".
    def create_graph_from_extensions(self, source_node, extGroup):
        '''Method to create or update the graph from the extension groups (e.g. reads overlapping with the source node that share the same extension)'''
        '''NB: add only the first read for each extension group, e.g. the read having the larger overlap'''
        self.add_node(source_node)
        for reads in extGroup.values():
            self.add_node(reads[0][0])
            self.add_edge((source_node, reads[0][0], len(source_node)-reads[0][1]))

    # Method "find_all_paths".
    def find_all_paths(self, start_node, end_node, path, all_paths):
        '''Method to return all the paths from the start_node to the end_node'''
        if start_node not in self.graph or end_node not in self.graph:
            return []
        # Path found from end_node to start_node.
        if end_node == start_node:
            path.reverse()
            all_paths.append(path)
            return
        # Traverse the graph to find the path from end_node to start_node.
        prev_nodes = []
        for node in self.graph:
            for neighbour in self.graph[node]:
                if end_node in neighbour:
                    prev_nodes.append(node)
        for node in prev_nodes:
            self.find_all_paths(start_node, node, path+[node], all_paths)
        return all_paths

    # Method "__repr__".
    def __repr__(self):
        return "Nodes graph: {}".format(self.graph)




#----------------------------------------------------
# getMostRepresentedKmer function
#----------------------------------------------------
def getMostRepresentedKmer(bamFile, region, kmerSize):
    """
    To get the most represented k-mer in a specific region of a BAM file. 

    Args:
        - bamFile: file
            BAM file of linked reads mapped onto the draft genome assembly
        - region: str
            region on the draft genome assembly from which to search for all alignments
            (format: "chr:posBeg-posEnd")
        - kmerSize: int
            k-mer size value

    Return/Output:
        - mostRepresentedKmer: str
            most represented k-mer's sequence
    """
    try:
        # Create a dictionary 'alignmentsOccurrencesDict' containing the number of occurrences of each alignment.
        alignmentsOccurrencesDict = {}

        # Region coordinates.
        chrID = region.split(':')[0]
        region_start = int(region.split(':')[1].split('-')[0])
        region_end = int(region.split('-')[-1])

        # Iterate over the alignments and update the 'alignmentsOccurrencesDict' dictionary.
        alignmentsFile = pysam.AlignmentFile(bamFile, "rb")
        for read in alignmentsFile.fetch(chrID, region_start, region_end):

            # Get the putative k-mer sequence.
            readSeq = str(read).split('\t')[9].split('array')[0]
            mappingPosition = int(str(read).split('\t')[3])
            putativeKmer = str(readSeq)[(region_start - mappingPosition +1):(region_start - mappingPosition + kmerSize +1)]

            # Update the 'alignmentsOccurrencesDict' dictionary.
            if putativeKmer in alignmentsOccurrencesDict:
                alignmentsOccurrencesDict[putativeKmer] += 1
            else:
                alignmentsOccurrencesDict[putativeKmer] = 1

        alignmentsFile.close()

        # Get the most represented k-mer (whose length is 'kmerSize').
        alignmentsOccurrencesFiltered = {k: v for k, v in alignmentsOccurrencesDict.items() if len(k) == kmerSize}
        if len(alignmentsOccurrencesFiltered) > 0:
            mostRepresentedKmer = max(alignmentsOccurrencesFiltered, key=alignmentsOccurrencesDict.get)
        else:
            mostRepresentedKmer = ""

        return mostRepresentedKmer

    except Exception as e:
        print("File 'helpers.py': Something wrong with the function 'getMostRepresentedKmer()'")
        print("Exception-{}".format(e))
        sys.exit(1)


#----------------------------------------------------
# updateGFAWithSolution function
#----------------------------------------------------
def updateGFAWithSolution(outDir, gfa_name, outputGFA, outputGFAFile):
    """
    To update the GFA, when a solution is found for a gap, with this solution.

    Args:
        - outDir: dir
            directory where the FASTA file containing all gap-filled sequences is located
        - gfa_name: str
            name of the GFA file
        - outputGFA: list
            list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the gap-filled sequence
        - outputGFAFile: file
            file of the output GFA (updated with the gap-filled sequence(s))

    Return/Output:
        - gapfillSeqFile: file
            file containing all the gap-filled sequences obtained with a good quality score [AB]
    """
    try:
        # Variables input.
        solutionName = outputGFA[0]
        seqLength = outputGFA[1]
        sequence = outputGFA[2]
        solution = outputGFA[3]
        pos_1 = outputGFA[4]
        pos_2 = outputGFA[5]
        leftName = solutionName.split(':')[0]
        rightName_w_orientation = solutionName.split(':')[1]
        rightName_wo_orientationSign = re.split('\+_|\-_', str(rightName_w_orientation))[0]
        rightOrientationSign = re.findall(r"\+_|\-_",str(rightName_w_orientation))[0].split('_')[0]
        rightName = rightName_wo_orientationSign + rightOrientationSign
        quality = outputGFA[6]

        print("Updating the GFA file with the solution: " + solutionName)

        # Save the gap-filled sequence to a file containing all gap-filled sequences.
        gapfillSeqFile = str(outDir) + "/" + str(gfa_name).split('.gfa')[0] + ".gapfilled_sequences.fasta"
        try:
            with open(gapfillSeqFile, "a") as seqFasta:
                seqFasta.write(">{} _ len.{}_qual.{} ".format(solutionName, seqLength, quality))
                seqFasta.write("\n" + sequence + "\n")
        except IOError as err:
            print("File 'helpers.py', function 'updateGFAWithSolution()': Unable to open or write to the file {}. \nIOError-{}".format(str(gapfillSeqFile), err))
            sys.exit(1)

        try:
            with open(outputGFAFile, "a") as f:
                # Add the gap-filled sequence (query seq) to GFA output ('Sequence' S line).
                out_gfa = gfapy.Gfa.from_file(outputGFAFile)
                out_gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(solutionName, seqLength, os.path.join(outDir, gapfillSeqFile)))

                # Write the two corresponding E lines ('Edges' lines) into GFA output.
                out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(leftName, solution, pos_1[0], pos_1[1], pos_1[2], pos_1[3]))
                out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(solution, rightName, pos_2[0], pos_2[1], pos_2[2], pos_2[3]))

                out_gfa.to_file(outputGFAFile)

        except IOError as err:
            print("File 'helpers.py', function 'updateGFAWithSolution()': Unable to open or write to the file {}. \nIOError-{}".format(str(outputGFAFile), err))
            sys.exit(1)

        return gapfillSeqFile
    
    except Exception as e:
        print("File 'helpers.py': Something wrong with the function 'updateGFAWithSolution()'")
        print("Exception-{}".format(e))
        sys.exit(1)

