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

"""Module 'helpers.py': Classes
#TODO
The module 'helpers.py' contains the classes and functions used in the script OLC.
"""

import os
import sys
import re
import subprocess
import gfapy
from gfapy.sequence import rc
from Bio import SeqIO
#from datetime import datetime      #TODO


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
# Graph class
#----------------------------------------------------
#TODO
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
# update_gfa_with_solution function
#----------------------------------------------------
# Function to update the GFA when a solution is found for a gap
def update_gfa_with_solution(outDir, gfa_name, output_for_gfa, gfa_output_file):
    """
    To update the GFA when a solution is found for a gap.

    Args:
        - outDir: dir
            directory where the FASTA file containing all gap-filled sequences is located
        - gfa_name: str
            name of the GFA file
        - output_for_gfa: list
            list containing the gap-filled sequence's name, as well as its length, its sequence, the number of solution found, the beginning and ending positions of the overlap and the quality of the sequence
        - gfa_output_file: file
            file of the output GFA (updated with the gap-filled sequence(s))

    Return/Output:
        - gapfillFile: file
            file containing all the gap-filled sequences found with a good quality score
    """
    try:
        # Variables input.
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

        # Save the gap-filled sequence to a file containing all gap-filled sequences.
        gapfillFile = gfa_name + ".gapfill_seq.fasta"
        with open(gapfillFile, "a") as seq_fasta:
            seq_fasta.write(">{} _ len_{}_qual_{} ".format(sol_name, length_seq, quality))
            seq_fasta.write("\n" + seq + "\n")

        with open(gfa_output_file, "a") as f:

            # Add the found seq (query seq) to GFA output (S line).
            out_gfa = gfapy.Gfa.from_file(gfa_output_file)
            out_gfa.add_line("S\t{}\t{}\t*\tUR:Z:{}".format(sol_name, length_seq, os.path.join(outDir, gapfillFile)))

            # Write the two corresponding E lines into GFA output.
            out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(s1, solution, pos_1[0], pos_1[1], pos_1[2], pos_1[3]))
            out_gfa.add_line("E\t*\t{}\t{}\t{}\t{}\t{}\t{}\t*".format(solution, s2, pos_2[0], pos_2[1], pos_2[2], pos_2[3]))

            out_gfa.to_file(gfa_output_file)

            return gapfillFile
    
    except Exception as e:
        print("\nFile 'helpers.py': Something wrong with the function 'update_gfa_with_solution()'")
        print("Exception-")
        print(e)
        sys.exit(1)

