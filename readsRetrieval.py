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

"""Module 'readsRetrieval.py': Reads Retrieval

The module 'readsRetrieval.py' enables to retrieve the reads whose barcode is observed in chunk regions surrounding the gap.
"""

from __future__ import print_function
import os
import re
import subprocess
import sys


#----------------------------------------------------
# retrieveReadsWithLRezQueryFastq function
#----------------------------------------------------
def retrieveReadsWithLRezQueryFastq(gfa_name, reads, index, allBarcodesLists, threads):
    """
    To retrieve the reads associated to the barcodes extracted on chunk regions, with `LRez query fastq`. 
    `LRez query fastq` enables to query a barcodes index and a fastq file to retrieve alignments containing the query barcodes.

    Args:
        - gfa_name: str
            name of the GFA file containing the gaps' coordinates
        - reads: file
            barcoded FASTQ file of linked reads
        - index: file
            index of barcodes
        - gapLabel: str
            label of the gap
        - allBarcodesLists: file
            file containing all the lists of the extracted barcodes (all the lists obtained for all gaps)
        - threads: int
            number of threads to use
    """
    try:
        # LRez query fastq.
        command = ["LRez", "query", "fastq", "--fastq", reads, "--index", index, "--collectionOfLists", allBarcodesLists, "--threads", str(threads), "--gzip"]
        queryFastqLog = str(gfa_name) + "_LRezQueryFastq.log"

        try:
            with open(queryFastqLog, "a") as log:
                subprocess.run(command, stderr=log)
        except IOError as err:
            print("File 'readsRetrieval.py', function 'retrieveReadsWithLRezQueryFastq()': Unable to open or write to the file {}. \nIOError-{}".format(str(queryFastqLog), err))
            sys.exit(1)

        # Remove the raw files obtained from `LRez query fastq`.
        if os.path.getsize(queryFastqLog) == 0:
            subprocess.run(["rm", queryFastqLog])

    except Exception as e:
        print("File 'readsRetrieval.py': Something wrong with the function 'retrieveReadsWithLRezQueryFastq()'")
        print("Exception-{}".format(e))
        sys.exit(1)

