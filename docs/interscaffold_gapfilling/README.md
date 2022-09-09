# Inter-scaffold gap-filling

Below is an example of the procedure used to perform inter-scaffold gap-filling with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.

## Input files

In the case of **inter-scaffold gap-filling**, you should have the following input files:
* A **matrix** file representing the links between scaffolds in tabular format
    * the matrix file contains the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds).
* A sequences file in **FASTA** format
    * the FASTA file contains scaffolds sequences.
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
# Matrix file: 'interscaffold_gapfilling.matrix'
# FASTA file: 'interscaffold_gapfilling.fasta'
# BAM file: 'interscaffold_gapfilling.bam'
# Indexed BAM file: 'interscaffold_gapfilling.bam.bai'
# FASTQ file (gzipped): 'interscaffold_gapfilling.fastq.gz'
```

## Generate the GFA file

To generate a GFA file from a matrix file containing the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds), you have to use the script **`matrix2gfa.py`** (located in utils/).
* The `matrix2gfa.py` script takes as input a file containing the number of common barcodes between all possibles pairs of scaffolds' extremities in tabular format (matrix) and converts it to a GFA file (GFA 2.0)
    * it is possible to specify the minimal number of links (minimal number of common barcodes) two scaffolds must share to consider a real link between them and so to add the corresponding gap to the GFA file: `-threshold` option
```
# Create a GFA file from the matrix file, with a minimal number of links between two scaffolds being of 10
../../utils/matrix2gfa.py -fa interscaffold_gapfilling.fasta -matrix interscaffold_gapfilling.matrix -out . -threshold 10
```
* interscaffold_gapfilling.fasta: FASTA file containing the sequences of the scaffolds
* interscaffold_gapfilling.matrix: Matrix file containing the number of common barcodes between all possibles pairs of scaffolds' extremities
* **Outputs**: 
    * GFA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa'
    * FASTA files containing the scaffolds sequences: 'interscaffold_gapfilling_[scaffoldID].scaffold.fasta'

### Build LRez barcode index

Prior to running MTG-Link, the LRez barcode index of the linked reads FASTQ file has to be built. This can be done with the following command:
```
# Build the LRez barcode index
LRez index fastq -f interscaffold_gapfilling.fastq.gz -o interscaffold_gapfilling.bci -g
```
* interscaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* **Output**:
    * LRez barcode index: 'interscaffold_gapfilling.bci'

### Run MTG-Link

MTG-Link can be run with the following command:  
```
# Run MTG-Link
../../mtglink.py DBG -gfa interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa -bam interscaffold_gapfilling.bam -fastq interscaffold_gapfilling.fastq.gz -index interscaffold_gapfilling.bci -t 4 -k 61 51 41 31 21
```
* interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa: GFA file obtained from `matrix2gfa.py`
* interscaffold_gapfilling.bam: input BAM file
* interscaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* interscaffold_gapfilling.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the outputs are located in the `MTG-Link_results/` directory
    * Output GFA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.
    * Optional output FASTA file: 'bad_solutions.fasta' (not output in this example)
        * another sequence file in FASTA format, that contains the assembled sequences found by MTG-Link but not returned in the output GFA file because of a bad quality score or a bad length (< 2*`-ext`).
    * `read_subsampling/` directory:
        * For each gap/target:
            * Barcodes file containing the barcodes observed in the gap/target flanking sequences ('.bxu')
            * Reads file in FASTQ format, containing the linked reads whose barcode is observed in the gap/target flanking sequences ('.bxu.fastq')
        * FOF barcodes file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa.flank10000.occ2.allBarcodesFiles.txt'
            * a file of file (FOF) containing a list of all the barcodes file. 
        * LOG file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa.flank10000.occ2.readSubsampling_summary.txt'
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
../../utils/gfa2tofasta.py -in ./MTG-Link_results/interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa -out .
```
* interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa: GFA file output by MTG-Link
* **Output**:
    * FASTA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink_assembly.fasta'

