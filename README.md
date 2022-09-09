# MTG-Link

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/mtglink)  
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


## What is MTG-Link ?

MTG-Link is a **local assembly** tool dedicated to **linked-reads**. It leverages barcode information from linked-reads to assemble **specific loci**.  
The locus of interest is defined by two coordinates on a draft genome, indicating its left and right flanking sequences, and the sequence in-between to be assembled is referred to as the target sequence.  
The main feature of MTG-Link is that it takes advantage of the linked-read barcode information to get a **subsample of reads** of interest for the **local assembly** of each sequence. It also automatically tests different parameters values and performs a **qualitative evaluation** of the obtained assemblies.  
MTG-Link can be used for various local assembly use cases, such as intra-scaffold and inter-scaffold gap-fillings, as well as alternative allele reconstruction of large insertion variants.

It takes as input a set of linked reads, the target flanking sequences and coordinates in GFA format (with the flanking sequences identified as ”segment” elements (S lines and the targets identified as ”gap” elements (G lines)) and an indexed BAM file obtained after mapping the linked reads onto the draft genome assembly.  
It outputs the set of assembled target sequences in Fasta format, as well as an assembly graph file in GFA format, complementing the input GFA file with the obtained target sequences.

Presently, it is directly compatible with the following linked-reads technologies, given the barcodes are reported using the BX:Z tag:
* 10X Genomics
* Haplotagging
* stLFR
* TELL-Seq

**MTG-Link** is a [Inria Genscale](https://team.inria.fr/genscale/) and [INRAE](https://www.inrae.fr/) / [INRIA](https://www.inria.fr/) tool developed by Anne Guichard.


## Installation

### External dependencies

* Biopython
* Gfapy
* Mummer
* Pathos
* Pysam
* LRez
* MindTheGap

For more information on the LRez and MindTheGap tools:
* LRez: <https://github.com/morispi/LRez>
* MindTheGap: <https://github.com/GATB/MindTheGap>

### Installation from conda

MTG-link is also distributed as a [Bioconda package](https://anaconda.org/bioconda/mtglink), which can be installed with:
```
conda install -c bioconda mtglink
```

### Installation from source

Clone the MTG-Link repository with:
```
git clone --recursive https://github.com/anne-gcd/MTG-Link.git
```

Please make sure you have installed all the external dependencies. Alternatively, you can install them via the "requirements.txt" file.  
To install a list of packages into a specified **conda** environment, do the following:  
`conda create --name <env> --file requirements.txt`  

### Testing the installation

You can test the installation of MTG-Link with the script `test.py` located in the `test/` directory, that runs MTG-Link on several test datasets, .   
Please run the following command to try out MTG-Link on these datasets:
```
cd test/
./test.py
```
To make sure MTG-Link is running properly, all tests should "Pass".  
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
* readsFile.fastq: Linked reads file. Warning: the barcode sequence must be in the header (BX:Z tag) TODO:below
* barcodeIndex.bci: File where to store the LRez barcode index

### 2) Run MTG-Link

MTG-Link can be run with the following command:  
```
mtglink.py DBG -gfa gfaFile.gfa -bam bamFile.bam -fastq readsFile.fastq.gz -index barcodeIndex.bci 
```
* gfaFile.gfa: GFA file containing the coordinates of the targets to fill
* bamFile.bam: BAM file of the linked reads mapped on the draft assembly. Warning: the associated .bai file must exist
* readsFile.fastq.gz: Linked reads file. Warning: the barcode sequence must be in the header (BX:Z tag). The FASTQ file must be gzipped
* barcodeIndex.bci: LRez barcode index of the FASTQ file

MTG-Link can be used for various local assembly use cases, such as intra-scaffold and inter-scaffold gap-fillings, as well as alternative allele reconstruction of large insertion variants.  
As users do not have the same input files depending on their use case, the `utils/` directory contains [scripts](./utils/README.md) to obtain the requested input GFA file from different input file formats.  
Besides, for each of these use cases, an example of the procedure to follow to perform local assembly with MTG-Link is detailed [here](./docs/UseCases.md).

#### Options

```
    -out OUTDIR             Output directory [default: ./MTG-Link_results]
    -line LINE              Line of GFA file input from which to start analysis (if not provided, start analysis from first line of GFA file input) [optional]
    -bxuDir BXUDIR          Directory where the FASTQ files containing the subsample of reads are located (1 file per target) (format of FASTQ files: xxx.bxu.fastq) [to provide if the read subsampling step has already been done for this dataset]
    -t THREADS              Number of threads to use for the Read Subsampling step [default: 1]
    -flank FLANKSIZE        Flanking sequences' size (bp) [default: 10000]
    -occ MINBARCOCC         Always remove barcodes for wich the number of occurrences in the union set from the two flanking sequences is smaller that this number [default: 2]
    -refDir REFDIR          Directory containing the reference sequences if any (1 file per target) [optional]
    -ext EXTSIZE            Size of the extension of the target on both sides (bp); determine start/end of local assembly [default: 500]
    -l MAXLENGTH            Maximum assembly length (bp) (it could be a bit bigger than the length of the target to fill OR it could be a very high length to prevent for searching indefinitely [default: 10000]
    -m MINLENGTH            Minimum assembly length (bp), by default 2*(-ext) bp [default: 1000]
    -k KMERSIZE             K-mer size(s) used for local assembly [default: [51, 41, 31, 21]]
    -a ABUNDANCETHRESHOLD   Minimal abundance threshold for solid k-mers [default: [3, 2]]
    --force                 To force search on all '-k' values provided
    -max-nodes MAXNODES     Maximum number of nodes in contig graph [default: 1000]
    -nb-cores NBCORES       Number of cores to use for the Local Assembly step (DBG assembly) [default: 1]
    -max-memory MAXMEMORY   Maximum memory for graph building (in MBytes) [default: 0]
    -verbose VERBOSITY      Verbosity level for DBG assembly [default: 0]
    --multiple              To return the assembled sequences even if multiple solutions are found (by default, if MTG-Link returns multiple solutions, we consider 'No Assembly' as it is not possible to know which one is the correct one)
```

**NB:** The `-refDir` directory should contain 1 file per target, e.g. 1 file per reference sequence. Besides, the files contained in `-refDir` should be formatted so that they contain the gap label. Thus, the prefix of these files and of their record ID should be for exemple:  
* Gap into scaffold: `<leftScaffoldName>_<coordLeftScaffoldStart>-<coordLeftScaffoldEnd>-L(+|-)_<rightScaffoldName>_<coordRightScaffoldStart>-<coordRightScaffoldEnd>-R(+|-)`
* Gap between scaffolds: `<ScaffoldName1>(+|-)_<scaffoldName2>(+|-)`

**NB:** When using the `--force` option, the `--multiple` option cannot be used, as otherwise it would filter unique solutions obtained with different `-k` values. 

<!--
## License
Please note that GATB-Core is distributed under Affero-GPL license.
-->

## Contact

To contact a developer, request help, or for any feedback on MTG-Link, please use the issue form of github: https://github.com/anne-gcd/MTG-Link/issues

You can see all issues concerning MTG-Link [here](https://github.com/anne-gcd/MTG-Link/issues).

If you do not have any github account, you can also send an email to Anne Guichard (<anne.guichard@irisa.fr>).

