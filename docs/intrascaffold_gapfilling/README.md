# Intra-scaffold gap-filling

Intra-scaffold gap-filling is defined as filling the gaps (regions of 'Ns') in scaffolds' sequences of a draft genome assembly.  
Below is an example of the procedure used to perform intra-scaffold gap-filling with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.


## How to generate the GFA file

**Input file**: A sequences file in **FASTA** format
* the FASTA file contains scaffolds sequences with 'Ns' regions (where 'Ns' regions will be treated as gaps/targets) (draft genome assembly).

**Scripts**: To generate a GFA file from a FASTA file containing 'Ns' regions, you have to use the scripts **`fasta2bed.py`** and **`bed2gfa.py`** (located in utils/).
* The `fasta2bed.py` script takes as input a FASTA file and converts it to a BED file containing the positions of 'Ns' for each scaffold
* The `bed2gfa.py` script converts a BED file containing the positions of 'Ns' for each scaffold to a GFA file (GFA 2.0) 
    * it is possible to filter the 'Ns' regions you want to treat as gaps/targets by:
        * their size (e.g. gap/target lengths): `-min` and `-max` options
        * the flanking contigs' sizes (for example, select only the 'Ns' regions whose flanking contigs' sizes are long enough to get enough barcodes): `-contigs` option
```
# Get a BED file containing the coordinates of the 'Ns' regions
../../utils/fasta2bed.py -fa intrascaffold_gapfilling.fasta -out .

# Create a GFA file with the gaps'/targets' coordinates of the 'Ns' regions of length 1000 bp and whose the flanking contigs' sizes are at least 10000 bp
../../utils/bed2gfa.py -bed intrascaffold_gapfilling_positions_Ns.bed -fa intrascaffold_gapfilling.fasta -out . -min 1000 -max 1000 -contigs 10000
```
* intrascaffold_gapfilling.fasta: FASTA file containing the sequences of the scaffolds (draft genome assemnly)
* intrascaffold_gapfilling_positions_Ns.bed: BED file obtained from `fasta2bed.py`
* **Outputs**: 
    * GFA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa'
    * FASTA files containing the left flanking region of the gap/target: 'intrascaffold_gapfilling_[scaffoldID]_[coordLeftFlankingSeq_start]-[coordLeftFlankingSeq_end].g1000.c10000.left.fasta'
    * FASTA files containing the right flanking region of the gap/target: 'intrascaffold_gapfilling_[scaffoldID]_[coordRightFlankingSeq_start]-[coordRightFlankingSeq_end].g1000.c10000.right.fasta'


## MTG-Link pipeline

```
# Build the LRez barcode index
LRez index fastq -f intrascaffold_gapfilling.fastq.gz -o intrascaffold_gapfilling.bci -g

# Run MTG-Link
../../mtglink.py DBG -gfa intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa -bam intrascaffold_gapfilling.bam -fastq intrascaffold_gapfilling.fastq.gz -index intrascaffold_gapfilling.bci -t 4 -k 61 51 41 31 21
```
* intrascaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa: GFA file obtained from `bed2gfa.py`
* intrascaffold_gapfilling.bam: input BAM file
* intrascaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* intrascaffold_gapfilling.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the outputs are located in the `MTG-Link_results/` directory
    * Output GFA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The output files and directory of MTG-Link are detailed [here](TODO)


## Get the updated FASTA file with assembled sequences

The GFA file output by MTG-Link is an assembly graph that represents the relationships between the input sequences and the target sequences (assemblies or gaps). The input FASTA file can be updated with the assembled sequences using this GFA file. The gaps of the GFA file are returned as 'Ns' regions in the output FASTA file.
```
# Get a BED file containing the coordinates of the 'Ns' regions
../../utils/gfa2tofasta.py -in ./MTG-Link_results/intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa -out .
```
* intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa: GFA file output by ***MTG-Link***
* **Output**:
    * FASTA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink_assembly.fasta'

