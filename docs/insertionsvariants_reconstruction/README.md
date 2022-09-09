# Alternative allele reconstruction of large insertion variants

Below is an example of the procedure used to perform alternative allele reconstruction of large insertion variants with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.

## Input files

In the case of **alternative allele reconstruction of large insertion variants**, you should have the following input files:
* A sequences file in **FASTA** format
    * the FASTA file contains scaffolds sequences (reference genome).
* The insertions coordinates in **VCF** format
    * the VCF file contains information on the large insertion variants (including the insertions coordinates).
* A set of linked reads in **FASTQ** format
    * the FASTQ file is a **barcoded** FASTQ file from linked reads. The barcode sequence must be in the header (BX:Z tag). For example, the *longranger basic* pipeline outputs a FASTQ file with barcode information attached.
    * the FASTQ file must be **gzipped**.
* An indexed **BAM** file obtained after mapping the linked reads onto the reference genome
    * the BAM file is a *Samtools* **indexed** BAM file. Each read in the BAM file has **barcode** information attached. For example, the *longranger* pipeline outputs an indexed BAM file containing position-sorted, aligned reads. Each read in this BAM file has Chromium barcode and phasing information attached.

**NB**: If you don't have the FASTQ file, you can generate it with the following command:
```
samtools bam2fq -T BX bamFile.bam > readsFile.fastq
gzip -c readsFile.fastq > readsFile.fastq.gz
```
* bamFile.bam: BAM file of the linked reads mapped on the reference genome. Warning: the associated .bai file must exist
* readsFile.fastq: Linked reads file, with the barcode sequence in the header (BX:Z tag)

### Example

```
# FASTA file: 'insertionsvariants_reconstruction.fasta'
# VCF file: 'insertionsvariants_reconstruction.vcf'
# BAM file: 'insertionsvariants_reconstruction.bam'
# Indexed BAM file: 'insertionsvariants_reconstruction.bam.bai'
# FASTQ file (gzipped): 'insertionsvariants_reconstruction.fastq.gz'
```

## Generate the GFA file

To generate a GFA file from a VCF file containing the insertions coordinates, you have to use the script **`vcf2gfa.py`** (located in utils/).
* The `vcf2gfa.py` converts a VCF file containing the insertions coordinates to a GFA file (GFA 2.0)
    * it is possible to filter the insertions you want to treat as targets (e.g. you want to reconstruct) by:
        * the flanking contigs' sizes (for example, select only the insertions whose flanking contigs' sizes are long enough to get enough barcodes): `-contigs` option
```
# Create a GFA file with the targets' coordinates of the insertions whose flanking contigs' sizes are at least 10000 bp
../../utils/vcf2gfa.py -vcf insertionsvariants_reconstruction.vcf -fa insertionsvariants_reconstruction.fasta -out . -contigs 10000
```
* insertionsvariants_reconstruction.vcf: VCF file containing the insertions coordinates
* insertionsvariants_reconstruction.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* **Outputs**: 
    * GFA file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa'
    * FASTA files containing the left flanking region of the target: 'insertionsvariants_reconstruction_[scaffoldID]_[coordLeftFlankingSeq_start]-[coordLeftFlankingSeq_end].g100.c10000.left.fasta'
    * FASTA files containing the right flanking region of the target: 'insertionsvariants_reconstruction_[scaffoldID]_[coordRightFlankingSeq_start]-[coordRightFlankingSeq_end].g100.c10000.right.fasta'

### Build LRez barcode index

Prior to running MTG-Link, the LRez barcode index of the linked reads FASTQ file has to be built. This can be done with the following command:
```
# Build the LRez barcode index
LRez index fastq -f insertionsvariants_reconstruction.fastq.gz -o insertionsvariants_reconstruction.bci -g
```
* insertionsvariants_reconstruction.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* **Output**:
    * LRez barcode index: 'insertionsvariants_reconstruction.bci'

### Run MTG-Link

MTG-Link can be run with the following command:  
```
# Run MTG-Link
../../mtglink.py DBG -gfa insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa -bam insertionsvariants_reconstruction.bam -fastq insertionsvariants_reconstruction.fastq.gz -index insertionsvariants_reconstruction.bci -t 4 -ext 100 -m 400 -k 51 41 31 21
```
* insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa: GFA file obtained from `vcf2gfa.py`
* insertionsvariants_reconstruction.bam: input BAM file
* insertionsvariants_reconstruction.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* insertionsvariants_reconstruction.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the outputs are located in the `MTG-Link_results/` directory
    * Output GFA file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.
    * Optional output FASTA file: 'bad_solutions.fasta' (not output in this example)
        * another sequence file in FASTA format, that contains the assembled sequences found by MTG-Link but not returned in the output GFA file because of a bad quality score or a bad length (< 2*`-ext`).
    * `read_subsampling/` directory:
        * For each target:
            * Barcodes file containing the barcodes observed in the target flanking sequences ('.bxu')
            * Reads file in FASTQ format, containing the linked reads whose barcode is observed in the target flanking sequences ('.bxu.fastq')
        * FOF barcodes file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa.flank10000.occ2.allBarcodesFiles.txt'
            * a file of file (FOF) containing a list of all the barcodes file. 
        * LOG file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa.flank10000.occ2.readSubsampling_summary.txt'
            * a tabular file with some information on the number of barcodes and reads extracted for eachtarget.
    * `local_assembly/` directory:
        * For each target:
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
        * For each target for which at least one assembled solution is found:
            * Sequence file in FASTA format, containing the target flanking sequences ('.contigs.fasta')

