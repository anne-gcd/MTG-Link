# MTG-Link

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/mtglink)  
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


## What is MTG-Link ?

MTG-Link is a novel **gap-filling** tool for draft genome assemblies, dedicated to **linked read** data.  
The main feature of MTG-Link is that it takes advantage of the linked-read barcode information to get a **subsample of reads** of interest for the **local assembly** of each gap. It also automatically tests different parameters values and performs a **qualitative evaluation** of the obtained solutions.

It takes as input a set of linked reads, a GFA file with gap coordinates and an indexed BAM file obtained after mapping the linked reads onto the draft assembly. It outputs the set of gap-filled sequences in FASTA format, as well as an assembly graph file in GFA format, containing the original contigs and the obtained gap-filled sequences of each gap, together with their overlapping relationships.

Presently, it is directly compatible with the following linked-reads technologies, given the barcodes are reported using the BX:Z tag:
* 10X Genomics
* Haplotagging
* stLFR
* TELL-Seq

**MTG-Link** is a [Inria Genscale](https://team.inria.fr/genscale/) and [INRAE](https://www.inrae.fr/) tool developed by Anne Guichard.


## Installation

### External dependencies

* Biopython
* Pathos
* Regex
* Gfapy
* Samtools
* Mummer
* Blast
* LRez
* MindTheGap

For more information on the LRez and MindTheGap tools:
* LRez: <https://github.com/morispi/LRez>
* MindTheGap: <https://github.com/GATB/MindTheGap>

### Installation from source

Clone the MTG-Link repository with:
```
git clone --recursive https://github.com/anne-gcd/MTG-Link.git
```

You can install the external dependencies via the **conda** package manager:  
`conda install -c bioconda samtools gfapy blast regex`  
`conda install -c conda-forge biopython pathos`   
`conda install -c bioconda/label/cf201901 mummer`  
`conda install -c bioconda mindthegap`  
`conda install -c bioconda lrez`  

Alternatively, you can install them via the requirements.txt file.  
To install a list of packages into a specified **conda** environment, do the following:  
`conda create --name <env> --file requirements.txt`  

For pysam, you have to install it using **pip**:  
`pip install pysam`  

### Installation from conda

MTG-link is also distributed as a [Bioconda package](https://anaconda.org/bioconda/mtglink), which can be installed with:
```
conda install -c bioconda mtglink
```

### Testing the installation

You can test the installation of MTG-Link with the script `test.sh`, that runs MTG-Link on several test datasets located in the `test/` directory.   
Please run the following command to try out MTG-Link on these datasets:
```
# Before running the tests, activate the conda environment containing the external dependencies or the MTG-Link bioconda package (called 'mtglink' in the example below)
conda activate mtglink
cd test/
./test.sh
```
In total, there are 10 tests. To make sure MTG-Link is running properly, all tests should "Pass".  
If at least one test "Fail", either MTG-Link is not installed properly or there is an unsolved issue in the code. In this case, please use the [issue form](https://github.com/anne-gcd/MTG-Link/issues) of GitHub.


## Running MTG-Link

### 1) Build LRez barcode index

Prior to running MTG-Link, the LRez barcode index of the linked reads FASTQ file has to be built. This can be done with one of the following command:
```
# The reads file is not gzipped.
LRez index fastq -f readsFile.fastq -o barcodeIndex.bci
# The reads file is gzipped.
LRez index fastq -f readsFile.fastq.gz -o barcodeIndex.bci -g
```
* readsFile.fastq: Linked reads file. Warning: the barcode sequence must be in the header (BX:Z tag)
* barcodeIndex.bci: File where to store the LRez barcode index

### 2) Run MTG-Link

The MTG-Link command line interface is composed of multiple parameters:
```
# The conda environment containing the external dependencies or the MTG-Link bioconda package is called 'mtglink'
conda activate mtglink
./mtglink.py --help

>>> usage: mtglink.py [-h] [-v] {DBG,IRO} ...

>>> Gapfilling with linked read data, using either a De Bruijn Graph (DBG) algorithm or an Iterative Read Overlap (IRO) algorithm

>>> positional arguments:
      {DBG,IRO}   MTGLink_module used for the Local Assembly step
        DBG       Gap-filling using a De Bruijn Graph (DBG) algorithm
        IRO       Gap-filling using an Iterative Read Overlap (IRO) algorithm

>>> optional arguments:
      -h, --help  show this help message and exit
      -v          show program's version number and exit
```

MTG-Link has two local assembly modules: **DBG** and **IRO**. For each module, the set of parameters is different. 

#### DBG module

MTG-Link can be run with the **DBG** module with the following command:  
```
./mtglink.py DBG -gfa gfaFile.gfa -bam bamFile.bam -fastq readsFile.fastq -index barcodeIndex.bci 
```
* gfaFile.gfa: GFA file containing the gaps to fill
* bamFile.bam: BAM file of the linked reads mapped on the draft assembly. Warning: the associated .bai file must exist
* readsFile.fastq: Linked reads file. Warning: the barcode sequence must be in the header (BX:Z tag)
* barcodeIndex.bci: LRez barcode index of the FASTQ file

##### Options

```
    -c CHUNK                Chunk size (bp) [default: 5000]
    -f FREQ                 Minimal frequence of barcodes observed in the union set from the two flanking gap sequences [default: 2]
    -out OUTDIR             Output directory [default: ./MTG-Link_results]
    -refDir REFDIR          Directory containing the reference sequences if any [optional]
    -line LINE              Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]
    -ext EXTENSION          Size of the extension of the gap on both sides (bp); determine start/end of gapfilling [default: 500]
    -l MAX_LENGTH           Maximum assembly length (bp) (it could be a bit bigger than the length of the gap to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]
    -k KMER_SIZE            K-mer size(s) used for gap-filling [default: [51, 41, 31, 21]]
    -a ABUNDANCE_THRESHOLD  Minimal abundance threshold for solid k-mers [default: [3, 2]]
    --force                 To force search on all '-k' values provided
    -t THREADS              Number of threads to use for the Read Subsampling step [default: 1]
    -max-nodes MAX_NODES    Maximum number of nodes in contig graph [default: 1000]
    -nb-cores NB_CORES      Number of cores to use for the Local Assembly step (DBG assembly) [default: 1]
    -max-memory MAX_MEMORY  Maximum memory for graph building (in MBytes) [default: 0]
    -verbose VERBOSITY      Verbosity level for DBG assembly [default: 0]
```

#### IRO module

MTG-Link can be run with the **IRO** module with the following command:  
```
./mtglink.py IRO -gfa gfaFile.gfa -bam bamFile.bam -fastq readsFile.fastq -index barcodeIndex.bci 
```
* gfaFile.gfa: GFA file containing the gaps to fill
* bamFile.bam: BAM file of the linked reads mapped on the draft assembly. Warning: the associated .bai file must exist
* readsFile.fastq: Linked reads file. Warning: the barcode sequence must be in the header (BX:Z tag)
* barcodeIndex.bci: LRez barcode index of the FASTQ file

##### Options

```
    -c CHUNK                Chunk size (bp) [default: 5000]
    -f FREQ                 Minimal frequence of barcodes observed in the union set from the two flanking gap sequences [default: 2]
    -out OUTDIR             Output directory [default: ./MTG-Link_results]
    -refDir REFDIR          Directory containing the reference sequences if any [optional]
    -line LINE              Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]
    -ext EXTENSION          Size of the extension of the gap on both sides (bp); determine start/end of gapfilling [default: 500]
    -l MAX_LENGTH           Maximum assembly length (bp) (it could be a bit bigger than the length of the gap to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]
    -s SEED_SIZE            Seed size used for indexing the reads (bp) [default: 10]
    -o MIN_OVERLAP          Minimum overlapping size (bp) [default: 20]
    -a ABUNDANCE_MIN        Minimal abundance(s) of reads used for gapfilling ; extension's groups having less than this number of reads are discarded from the graph [default: [3, 2]]
    -dmax MAX_SCORE         Maximum number of gaps/substitutions allowed in the inexact overlap between reads [default: 2]
    -t THREADS              Number of threads to use for the Read Subsampling step [default: 1]
```


## User Manual

MTG-Link takes as input a GFA file (GFA 2.0) with gap coordinates, a set of linked reads, an indexed BAM file obtained after mapping the linked reads onto the draft assembly and an indexed FASTQ file.  
It outputs the results in a GFA file (GFA 2.0), containing the original contigs and the obtained gap-filled sequences of each gap, together with their overlapping relationships. It also returns the set of gap-filled sequences in a FASTA file.

The qualitative evaluation of MTG-Link can be performed either with the corresponding reference sequences of the gaps (`-refDir`) or with the flanking contigs sequences.  
**NB:** The files contained in `-refDir` should be formatted so that they contain the gap label. Thus, the prefix of these files and of their record ID should be for ex:  
* Gap into scaffold: `<scaffoldName>_<coordFlankingScaffoldStart>-<coordFlankingScaffoldEnd>-L(+|-)_<scaffoldName>_<coordFlankingScaffoldStart>-<coordFlankingScaffoldEnd>-R(+|-)`
* Gap between scaffolds: `<scaffoldName1>(+|-)_<scaffoldName2>(+|-)`
#TODO: Move to "Preparing input files"

### Preparing input files

#### GFA file

The **GFA** (Graphical Fragment Assembly) file is a *tab-delimited* file containing the gap coordinates. The expected format is a [GFA 2.0](http://gfa-spec.github.io/GFA-spec/GFA2.html) format:  
```
<header>   <- H {VN:Z:2.0}
<segment>  <- S <sid1> <slen1> * UR:Z:<path_to_fasta_sequence>
<segment>  <- S <sid2> <slen2> * UR:Z:<path_to_fasta_sequence>
<gap>      <- G (* | <gid>) <sid1(+|-)> <sid2(+|-)> <dist> (* | <var>)
```
* \<sid\>: ID of the gap flanking sequence
* \<slen\>: Length of the gap flanking sequence (bp)
* \<path_to_fasta_sequence\>: Path of the FASTA file containing the sequence of the gap flanking sequence
* \<gid\>: ID of the gap
* \<dist\>: Estimated gap length (bp)
* \<var\>: Variance of the gap length estimate

The format of the `<sid>` from the **S-lines** depends on whether the gap is located into the scaffold or between scaffolds:  
* Gap into scaffold: `<scaffoldName>_<coordFlankingScaffoldStart>-<coordFlankingScaffoldEnd>-(L|R)`  
* Gap between scaffolds: `<scaffoldName>`

The **FASTA** files containing the segment sequences (indicated after `UR:Z:` in the GFA file) should have record ID that starts with the ID of the corresponding segment (`<sid>`).  

The **G-lines** describe a gap edge, that gives the estimated gap distance between the two segment sequences and the variance of that estimate.  
The gap is between the first segment at left `<sid1(+|-)>` and the second segment at right `<sid2(+|-)>` where the segments are oriented according to their sign indicators `(+|-)`.  
The `<dist>` field gives the expected distance between the first and second segment in their respective orientations, or 0 is this expected distance is unknown.  

Example of a GFA file containing a gap into the scaffold 8 (gap located between the positions 3152098 bp and 3153098 bp ; estimated length of 1000 bp):
```
    H	VN:Z:2.0
    S	8_0-3152098-L	3152098	*	UR:Z:8_0-3152098.g1000.c5000.left.fasta
    S	8_3153098-6305195-R	3152097	*	UR:Z:8_3153098-6305195.g1000.c5000.right.fasta
    G	*	8_0-3152098-L+	8_3153098-6305195-R+	1000	*
```

How to obtain a GFA file:  
* If you have a file containing the paths between scaffolds, you can use the **`paths2gfa.py`** script (in the `utils/` directory).  
  Format of a path: `<int:nb_scaffolds>****<sid1(f|r)>+<sid2(f|r)>`
* If you have a FASTA file with sequences containing 'Ns' regions (where 'Ns' regions will be treated as gaps), you can use the **`fasta2bed.py`** and **`bed2gfa.py`** scripts (in the `utils/` directory).
* If you have a file containing the links between the ends of the scaffolds in tabular format (e.g. a matrix), you can use the **`matrix2gfa.py`** script (in the `utils/` directory).

#### BAM file

The **BAM** file is a *Samtools* **indexed** BAM file, obtained after mapping the linked reads onto the assembly. Each read in the BAM file has **barcode** information attached.     
For example, the *longranger* pipeline outputs an indexed BAM file containing position-sorted, aligned reads. Each read in this BAM file has Chromium barcode and phasing information attached.  

#### FASTQ file

The **FASTQ** file is a **barcoded** FASTQ file from linked reads. The barcode sequence must be in the header (BX:Z tag).  
For example, the *longranger basic* pipeline outputs a FASTQ file with barcode information attached.  
**ATTN:** The FASTQ file must be gzipped. 

#### Index file

The **index** file contains the barcode index of the FASTQ file.  
To get the index file, you need to use the subcommand `LRez index fastq` of the tool [LRez](https://github.com/morispi/LRez):

### Output files

MTG-Link outputs several files:  
* an assembly graph file (`_mtglink.gfa`) in GFA format. It contains the original contigs and the obtained gap-filled sequences of each gap, together with their overlapping relationships. 
* a sequence file (`.gapfilled_sequences.fasta`) in FASTA format. It contains the set of gap-filled sequences.
* [optional] another sequence file ('`bad_solutions.fasta`') in FASTA format. It contains the inserted sequences found by MTG-Link but not returned in the output GFA file because of a bad quality score or a bad length (< 2*`-ext`).

There is also a `read_subsampling/` directory, with:  
* a barcodes file (`.bxu`), containing the barcodes observed in the gap flanking sequences.
* a FOF (`.allBarcodesFiles.txt`)), containing a list of all the barcodes file. 
* a reads file (`.bxu.fastq`) in FASTQ format. It contains the linked reads whose barcode is observed in the gap flanking sequences.
* a log file (`.readSubsampling_summary.txt`), a tabular file with some information on the number of barcodes and reads extracted for each gap.

There is also a `local_assembly/` directory, with:  
* For the **DBG module**:  
    * a breakpoint file (`.bkpt.fasta`) in FASTA format. It contains the breakpoint sequences (e.g the start and stop k-mers) used for gap-filling.
    * a log file (`.info.txt`), a tabular file with some information about the filling process for each breakpoint.
    * a sequence file (`.insertions_quality.fasta`) in FASTA format. It contains the inserted/gap-filled sequences that were successfully assembled, with their qualitative scores. 
* For the **IRO module**:  
    * a breakpoint file (`.kmersStartStop.fa`) in FASTA format. It contains the breakpoint sequences (e.g the start and stop k-mers) used for gap-filling.
    * a log file (`_IROresults.log`), a text file with extended (but not completed) sequences and some information on why the gap-filling was not completed for these sequences.
    * a sequence file (`.insertions_quality.fasta`) in FASTA format. It contains the inserted/gap-filled sequences that were successfully assembled, with their qualitative scores. 

There is also a `qual_evaluation/` directory, with:  
* a log file (`.ref_qry.log`), a text file with some information on the input files used for the alignment. 
* an alignment file (`.nucmerAlignments.stats`), a tabular file with some information on the alignment performed between the gap-filled sequences and a reference sequence.

There is also a `contigs/` directory, when the qualitative evaluation is performed with the flanking contigs, with:  
* a sequence file (`.contigs.fasta`) in FASTA format. It contains the gap flanking sequences.


## Description of MTG-Link

For each gap, MTG-Link relies on a **three-step pipeline** to gap-fill the gap:  
* The first step uses the barcode information of the linked read dataset to get a **subsample of reads** of potential interest for gap-filling.  
* The second step performs **local assembly** using this subsample of linked reads. Two different assembly algorithms are implemented and can be interchangeably used. The first one, called hereafter the **De Bruijn Graph (DBG) algorithm**, uses a de Bruijn graph data structure, and the second one, called the **Iterative Read Overlap (IRO) algorithm**, is based on on-the-fly computations of read overlaps.  
* The third step evaluates the obtained gap-filled sequence and annotates it with a **quality score**.  

In order to speed up the process, MTG-Link uses a trivial **parallelization** scheme by giving each gap to a separate thread.   

![MTG-Link_pipeline](doc/images/Overview_MTGLink_gapfilling_pipeline.png)

### a) Read subsampling

For each gap, MTG-Link extracts the linked reads whose barcode is observed in the chunk regions surrounding the gap, using the thirdparty tool **LRez**.  
The chunk region size can be defined by the user, the default value being 5,000 bp.  
To increase specificity, we keep only the barcodes for which the number of occurrences in the union set from the two flanking sequences is larger than a user-defined parameter -f (by default 2).  

The linked reads extracted during this step constitute the read subsample used in the local assembly step.

### b) Local assembly

To fill the gap between two contigs, a local assembly is performed using the subsample of linked reads obtained during the "Read subsampling" step.  
To be able to further evaluate the obtained gap-filled sequence, we extend the gap on both sides by `-ext` bp (by default 500 bp). Thus, MTG-Link will perform the local assembly between the sequences surrounding the extended gap, e.g. from the k-mer START (source) to the k-mer STOP (target).   

Two assembly algorithms can be used during this step: the **DBG algorithm** or the **IRO algorithm**.

#### Module DBG

The DBG algorithm is performed with the ***fill*** module of the software **MindTheGap**. The ***fill*** module of MindTheGap is an efficient local assembly module that relies on a **De Bruijn graph** data structure to represent the input read sequences. MindTheGap will try to find a path in the de Bruijn graph from the k-mer START (source) to the k-mer STOP (target). The gap-filling is performed in both forward and reverse orientation.  
In this module, two parameters have major impacts on the quality of the assembly: 
* The **k-mer size** (**`-k`**)  
* The **k-mer abundance threshold** (**`-a`**) for including a k-mer in the graph (solid k-mer threshold)  

These parameters are usually set in accordance with the expected sequencing depth. In the case of MTG-link, the latter may vary depending on the efficiency of the barcode-based subsampling step. Hence for higher sensitivity, MTG-Link automatically tests different values for these two parameters, starting with the highest ones and decreasing the values if no inserted sequence is found.

#### Module IRO

The IRO algorithm is based on **on-the-fly computations of read overlaps and iterative extensions** of the current assembly sequence. Overlapping reads are reads whose prefix (or reverse complement of the suffix) aligns with the suffix of the current assembly sequence with at most **`-dmax`** differences (including
substitutions and indels) over at least **`-o`** bp. These overlaps are found using a **seed-and-extend** schema, combining a **seed indexing** with a hash table and a **banded dynamic programming semi-global alignment algorithm**.  
At each iteration, several possible extensions may be found, due to sequencing errors and/or repeats. To avoid including sequencing errors, only extensions that are supported by a minimum number of reads (parameter **`-a`**, by default 2) are considered. Then, another extension phase begins.  
In the following cases, the algorithm **backtracks** and tries other extensions previously encountered but not yet explored:  
* When no overlapping read is found  
* If there is no extension shared by a sufficient number of reads  
* If the maximal assembled sequence size (user defined parameter) is reached  

Finally, if during an extension phase, the k-mer STOP is found, the assembly sequence is returned and the exploration ends.

### c) Qualitative evaluation

MTG-Link assigns a **quality score** to each gap-filled sequence obtained during the "Local assembly" step, that might help filtering out putative erroneous sequences:
* If a **reference sequence** is provided (**`-refDir`**):  1-letter score X<sub>1</sub> with X = [A, B, C, D]
    * X<sub>1</sub>: alignment to the entire reference sequence
* Using the **gap flanking contigs** information (**`-ext`**): 2-letters score X<sub>1</sub>X<sub>2</sub> with X = [A, B, C, D]
    * X<sub>1</sub>: alignment to the left flanking sequence (left 'ext' sequence)
    * X<sub>2</sub>: alignment to the right flanking sequence (right 'ext' sequence)

MTG-Link selects the gap-filled sequences with a score **[AB]** (reference sequence provided) or with a score **[AB]{2}** (using the flanking contigs information.  
To have a good quality score:  
* When comparing to the flanking contigs:  
    * The gap-filled sequence must be larger than twice `-ext` bp  
    * The gap-filled sequence must align on at least 90% of the lengths of the gap flanking sequences  
* When comparing to the flanking contigs:  
    * The length of the gap-filled sequence must be +-10% of the reference length  
    * The gap-filled sequence must align on at least 90% of the reference length

### Output

MTG-Link returns the results in a **GFA** file (GFA 2.0), containing the original contigs and the obtained gap-filled sequences of each gap, together with their overlapping relationships. It also returns the set of gap-filled sequences in a **FASTA** file. 

<!--
## License
Please note that GATB-Core is distributed under Affero-GPL license.
-->

## Contact

To contact a developer, request help, or for any feedback on MTG-Link, please use the issue form of github: https://github.com/anne-gcd/MTG-Link/issues

You can see all issues concerning MTG-Link [here](https://github.com/anne-gcd/MTG-Link/issues).

If you do not have any github account, you can also send an email to Anne Guichard (<anne.guichard@irisa.fr>).

