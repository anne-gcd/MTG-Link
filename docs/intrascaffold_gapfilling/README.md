# Intra-scaffold gap-filling

Below is an example of the procedure used to perform intra-scaffold gap-filling with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.

## Input files

In the case of **intra-scaffold gap-filling**, you should have the following input files:
* A draft genome assembly in **FASTA** format
    * the FASTA file contains scaffolds sequences with 'Ns' regions (where 'Ns' regions will be treated as gaps/targets).
* A set of linked reads in **FASTQ** format
    * the FASTQ file is a **barcoded** FASTQ file from linked reads. The barcode sequence must be in the header (BX:Z tag). For example, the *longranger basic* pipeline outputs a FASTQ file with barcode information attached.
    * the FASTQ file must be **gzipped**.
* An indexed **BAM** file obtained after mapping the linked reads onto the draft genome assembly
    * the BAM file is a *Samtools* **indexed** BAM file. Each read in the BAM file has **barcode** information attached. For example, the *longranger* pipeline outputs an indexed BAM file containing position-sorted, aligned reads. Each read in this BAM file has Chromium barcode and phasing information attached.

**NB**: If you don't have the FASTQ file, you can generate it with the following command:
```
samtools bam2fq -T BX bamFile.bam > readsFile.fastq
gzip -c readsFile.fastq > readsFile.fastq.gz
```
* bamFile.bam: BAM file of the linked reads mapped on the draft assembly. Warning: the associated .bai file must exist
* readsFile.fastq: Linked reads file, with the barcode sequence in the header (BX:Z tag)

### Example

```
# FASTA file: 'intrascaffold_gapfilling.fasta'
# BAM file: 'intrascaffold_gapfilling.bam'
# Indexed BAM file: 'intrascaffold_gapfilling.bam.bai'
# FASTQ file (gzipped): 'intrascaffold_gapfilling.fastq.gz'
```

## Generate the GFA file

To generate a GFA file from a FASTA file containing 'Ns' regions, you have to use the scripts **`fasta2bed.py`** and **`bed2gfa.py`** (located in utils/).
* The `fasta2bed.py` script takes as input a FASTA file and converts it to a BED file containing the positions of 'Ns' for each scaffold
* The `bed2gfa.py` script converts a BED file containing the positions of 'Ns' for each scaffold to a GFA file (GFA 2.0) 
    * it is possible to filter the 'Ns' regions you want to treat as gaps/targets by:
        * their size (e.g. gap/target lengths): `-min` and `-max` options
        * the flanking contigs' sizes (for example, select only the 'Ns' regions whose flanking contigs' sizes are long enough to get enough barcodes): `-contigs` option
```
# Get a BED file containing the coordinates of the 'Ns' regions
../../utils/fasta2bed.py -fa intrascaffold_gapfilling.fasta -out .
```
* intrascaffold_gapfilling.fasta: FASTA file containing the sequences of the scaffolds obtained from the draft assembly
* **Output**:
    * BED file: 'intrascaffold_gapfilling_positions_Ns.bed'
```
# Create a GFA file with the gaps'/targets' coordinates of the 'Ns' regions of length 1000 bp and whose the flanking contigs' sizes are at least 10000 bp
../../utils/bed2gfa.py -bed intrascaffold_gapfilling_positions_Ns.bed -fa intrascaffold_gapfilling.fasta -out . -min 1000 -max 1000 -contigs 10000
```
* intrascaffold_gapfilling_positions_Ns.bed: BED file obtained from `fasta2bed.py`
* intrascaffold_gapfilling.fasta: FASTA file containing the sequences of the scaffolds obtained from the draft assembly
* **Outputs**: 
    * GFA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa'
    * FASTA files containing the left flanking region of the gap/target: 'intrascaffold_gapfilling_[scaffoldID]_[coordLeftFlankingSeq_start]-[coordLeftFlankingSeq_end].g1000.c10000.left.fasta'
    * FASTA files containing the right flanking region of the gap/target: 'intrascaffold_gapfilling_[scaffoldID]_[coordRightFlankingSeq_start]-[coordRightFlankingSeq_end].g1000.c10000.right.fasta'

### Build LRez barcode index

Prior to running MTG-Link, the LRez barcode index of the linked reads FASTQ file has to be built. This can be done with the following command:
```
# Build the LRez barcode index
LRez index fastq -f intrascaffold_gapfilling.fastq.gz -o intrascaffold_gapfilling.bci -g
```
* intrascaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* **Output**:
    * LRez barcode index: 'intrascaffold_gapfilling.bci'

### Run MTG-Link

MTG-Link can be run with the following command:  
```
# Run MTG-Link
../../mtglink.py DBG -gfa intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa -bam intrascaffold_gapfilling.bam -fastq intrascaffold_gapfilling.fastq.gz -index intrascaffold_gapfilling.bci -t 4 -k 61 51 41 31 21
```
* intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa: GFA file obtained from `bed2gfa.py`
* intrascaffold_gapfilling.bam: input BAM file
* intrascaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* intrascaffold_gapfilling.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the outputs are located in the `MTG-Link_results/` directory
    * Output GFA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.
    * Optional output FASTA file: 'bad_solutions.fasta' (not output in this example)
        * another sequence file in FASTA format, that contains the assembled sequences found by MTG-Link but not returned in the output GFA file because of a bad quality score or a bad length (< 2*`-ext`).
    * `read_subsampling/` directory:
        * For each gap/target:
            * Barcodes file containing the barcodes observed in the gap/target flanking sequences ('.bxu')
            * Reads file in FASTQ format, containing the linked reads whose barcode is observed in the gap/target flanking sequences ('.bxu.fastq')
        * FOF barcodes file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa.flank10000.occ2.allBarcodesFiles.txt'
            * a file of file (FOF) containing a list of all the barcodes file. 
        * LOG file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.gfa.flank10000.occ2.readSubsampling_summary.txt'
            * a tabular file with some information on the number of barcodes and reads extracted for each gap/target.
    * `local_assembly/` directory:
        * For each gap/target:
            * Breakpoint file in FASTA format, containing the breakpoint sequences (e.g the start and stop k-mers) used for the local assembly ('.bkpt.fasta')
            * Raw sequence file in FASTA format, containing the assembled sequences found by the local assembly step of MTG-Link ('.insertions.fasta')
            * Filtered sequence file in FASTA format, containing the assembled sequences found by the local assembly step of MTG-Link and having an assembly length larger than the minimum assembly length required ('.insertions_filtered.fasta')
            * Qualitatively evaluated sequence file in FASTA format, containing the filtered assembled sequences found by the local assembly step of MTG-Link with their qualitative evaluation score ('.insertions_filtered_quality.fasta')
            * LOG file: a tabular file with some information about the filling process for each breakpoint ('.info.tx')
    * `qual_evaluation/` directory:
        * For each assembled target:
            * LOG file: a text file with some information on the input files used for the alignment ('.ref_qry.log')
            * Alignment file: a tabular file with some metrics on the alignment performed between the assembled sequences and (a) reference sequence(s) ('.nucmerAlignments.stats')
                * and the unsorted alignment file ('.nucmerAlignments.stats.unsorted')
    * `contigs/` directory: [when the qualitative evaluation is performed with the flanking contigs]
        * For each gap/target for which at least one assembled solution is found:
            * Sequence file in FASTA format, containing the gap/target flanking sequences ('.contigs.fasta')

### Get the updated FASTA file with assembled sequences

The GFA file output by MTG-Link is an assembly graph that represents the relationships between the input sequences and the target sequences (assemblies or gaps). The input FASTA file can be updated with the assembled sequences using this GFA file. The gaps of the GFA file are returned as 'Ns' regions in the output FASTA file.
```
# Get a BED file containing the coordinates of the 'Ns' regions
../../utils/gfa2tofasta.py -in ./MTG-Link_results/intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa -out .
```
* intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa: GFA file output by MTG-Link
* **Output**:
    * FASTA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink_assembly.fasta'

