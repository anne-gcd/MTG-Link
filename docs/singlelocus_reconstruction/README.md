# A single locus of interest

We are interested here in a particular locus on a reference genome. We have re-sequenced an individual from the same species with a linked-read technology and we want to assemble this specific locus from the sequencing reads of this individual. 

The locus of interest is defined by two coordinates on a chromosome or scaffold of the reference genome, therefore we dispose of a BED file (locus coordinates) and a sequence file (reference genome) to generate the input GFA file for MTG-Link.


## How to generate the GFA file

**Input file**: A sequences file in **FASTA** format
* the FASTA file contains the chromosome or scaffold sequence (reference genome).

**Scripts**: To generate a GFA file from a BED file containing the locus coordinates, you have to use the script **`bed2gfa.py`** (located in `utils/`).
* The `bed2gfa.py` script converts a BED file containing the locus coordinates to a GFA file (GFA 2.0) 

```
# Generate the input GFA file
../../utils/bed2gfa.py -bed myLocus.bed -fa referenceGenome.fasta -out myLocus.gfa
```
* referenceGenome.fasta: FASTA file containing the chromosome or scaffold sequence (reference genome)
* myLocus.bed: BED file containing the locus coordinates
* **Outputs**: 
    * GFA file: 'myLocus.gfa'
    * FASTA files containing the left flanking region of the gap/target: 'referenceGenome_[scaffoldID]_[coordLeftFlankingSeq_start]-[coordLeftFlankingSeq_end].g[locusLength].left.fasta'
    * FASTA files containing the right flanking region of the gap/target: 'referenceGenome_[scaffoldID]_[coordRightFlankingSeq_start]-[coordRightFlankingSeq_end].g[locusLength].right.fasta'


## MTG-Link pipeline

```
# Build the LRez barcode index (if not already done)
LRez index fastq -f readsFile.fastq.gz -o barcodeIndex.bci -g

# Run MTG-Link
../../mtglink.py DBG -gfa myLocus.gfa -bam bamFile.bam -fastq readsFile.fastq.gz -index barcodeIndex.bci -t 4 -k 61 51 41 31 21
```
* readsFile.fastq.gz: Linked-reads file (must be gzipped), with the barcode sequence in the header (BX:Z tag)
* myLocus.gfa: GFA file obtained from `bed2gfa.py`
* bamFile.bam: input BAM file
* barcodeIndex.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the main outputs of MTG-Link are the following:
    * Output GFA file: 'myLocus_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'myLocus.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The outputs of MTG-Link are detailed [here](../input-output_files.md)


## Get the updated FASTA file with assembled sequences

The GFA file output by MTG-Link is an assembly graph that represents the relationships between the locus flanking sequences and the locus sequence (assemblies or gaps). The input FASTA file can be updated with the assembled locus sequence using this GFA file. The gaps of the GFA file are returned as 'N's regions in the output FASTA file.
```
# Update the input FASTA file with assembled sequences
../../utils/gfa2tofasta.py -in ./MTG-Link_results/myLocus_mtglink.gfa -out .
```
* myLocus_mtglink.gfa: GFA file output by ***MTG-Link***
* **Output**:
    * FASTA file: 'myLocus_mtglink_assembly.fasta'

