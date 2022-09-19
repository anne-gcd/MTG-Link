# Alternative allele reconstruction of large insertion variants

Alternative allele reconstruction of large insertion variants is defined as reconstructing the sequence of the insertion variant present in the alternative allele.  
Below is an example of the procedure used to perform alternative allele reconstruction of large insertion variants with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.


## How to generate the GFA file

**Input file**: The insertions coordinates in **VCF** format
* the VCF file contains information on the large insertion variants (including the insertions coordinates).

**Scripts**: To generate a GFA file from a VCF file containing the insertions coordinates, you have to use the script **`vcf2gfa.py`** (located in utils/).
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


## MTG-Link pipeline

```
# Build the LRez barcode index
LRez index fastq -f insertionsvariants_reconstruction.fastq.gz -o insertionsvariants_reconstruction.bci -g

# Run MTG-Link
../../mtglink.py DBG -gfa insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa -bam insertionsvariants_reconstruction.bam -fastq insertionsvariants_reconstruction.fastq.gz -index insertionsvariants_reconstruction.bci -t 4 -ext 100 -m 400 -k 51 41 31 21
```
* insertionsvariants_reconstruction.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.gfa: GFA file obtained from `vcf2gfa.py`
* insertionsvariants_reconstruction.bam: input BAM file
* insertionsvariants_reconstruction.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the outputs are located in the `MTG-Link_results/` directory
    * Output GFA file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'insertionsvariants_reconstruction_insertions_extension_50_contigs_10000.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The outputs of MTG-Link are detailed [here](../input-output_files.md)

