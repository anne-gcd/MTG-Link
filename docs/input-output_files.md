# Input files

For each use case, the following files are common input files required by MTG-Link:
* A sequences file in **FASTA** format
    * the FASTA file contains scaffolds sequences.
* A set of linked reads in **FASTQ** format
    * the FASTQ file is a **barcoded** FASTQ file from linked reads. The barcode sequence must be in the header (BX:Z tag). For example, the *longranger basic* pipeline outputs a FASTQ file with barcode information attached.
    * the FASTQ file must be **gzipped**.
* An indexed **BAM** file obtained after mapping the linked reads onto the reference genome
    * the BAM file is a *Samtools* **indexed** BAM file. Each read in the BAM file has **barcode** information attached. For example, the *longranger* pipeline outputs an indexed BAM file containing position-sorted, aligned reads. Each read in this BAM file has Chromium barcode and phasing information attached.

**NB**: If you don't have the FASTQ file, you can generate it with the following command:
```
samtools bam2fq -T BX bamFile.bam | gzip > readsFile.fastq.gz
```
* bamFile.bam: BAM file of the linked reads mapped on the reference genome. Warning: the associated .bai file must exist
* readsFile.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)

Besides, depending on the use case, an additional file may be required and is specific to each use case (see documentation for each use case).


# Output files and directories

The main outputs of MTG-Link are the following:
* Output GFA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000_mtglink.gfa'
    * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
* Output FASTA file: 'intrascaffold_gapfilling_gaps_1000-1000_contigs_10000.assembled_sequences.fasta'
    * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

In the case of an assembly not validated by the qualitative evaluation step for at least one target, MTG-Link also outputs the following file:
* Optional output FASTA file: 'bad_solutions.fasta' (not output in this example)
    * another sequence file in FASTA format, that contains the assembled sequences found by MTG-Link but not returned in the output GFA file because of a bad quality score or a bad length (< 2*`-ext`).

Besides, MTG-Link outputs directories containing intermediate files obtained from each step of the pipeline:
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

