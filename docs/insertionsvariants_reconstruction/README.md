# Alternative allele reconstruction of large insertion variants

Alternative allele reconstruction of large insertion variants is defined as reconstructing the insertions sequences of all insertions calls from a VCF file.  
Below is an example of the procedure used to perform alternative allele reconstruction of large insertion variants with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.


## How to generate the GFA file

**Input file**: 
* The insertions coordinates in **VCF** format
    * the VCF file contains information on the large insertion variants (including the insertions coordinates).

**Scripts**: To generate a GFA file from a VCF file containing the insertions coordinates, you have to use the script **`vcf2gfa.py`** (located in `utils/`).
* The `vcf2gfa.py` converts a VCF file containing the insertions coordinates to a GFA file (GFA 2.0)
    * only the insertions calls defined with the 'INS' field are treated
    * it is possible to filter the insertions you want to treat as targets (e.g. you want to reconstruct) by:
        * the flanking contigs' sizes (for example, select only the insertions whose flanking contigs' sizes are long enough to get enough barcodes): `-contigs` option
```
# Create a GFA file with the targets' coordinates of the insertions whose flanking contigs' sizes are at least 10000 bp
../../utils/vcf2gfa.py -vcf vcfFile.vcf -fa referenceGenome.fasta -out insertionSequenceCoordinates.gfa -contigs 10000
```
* `vcfFile.vcf`: VCF file containing the insertions coordinates
* `referenceGenome.fasta`: FASTA file of the reference genome
* **Outputs**: 
    * GFA file: `insertionSequenceCoordinates.gfa`
    * FASTA files containing the left flanking region of the target: `referenceGenome_[scaffoldID]_[coordLeftFlankingSeq_start]-[coordLeftFlankingSeq_end].g100.left.fasta`
    * FASTA files containing the right flanking region of the target: `referenceGenome_[scaffoldID]_[coordRightFlankingSeq_start]-[coordRightFlankingSeq_end].g100.right.fasta`


## MTG-Link pipeline

```
# Build the LRez barcode index
LRez index fastq -f readsFile.fastq.gz -o barcodeIndex.bci -g

# Run MTG-Link
../../mtglink.py DBG -gfa insertionSequenceCoordinates.gfa -bam bamFile.bam -fastq readsFile.fastq.gz -index barcodeIndex.bci -ext 100 -m 400 -k 51 41 31 21
```
* `readsFile.fastq.gz`: Linked-reads file (must be gzipped), with the barcode sequence in the header (BX:Z tag)
* `insertionSequenceCoordinates.gfa`: GFA file obtained from `vcf2gfa.py`
* `bamFile.bam`: input BAM file
* `barcodeIndex.bci`: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the main outputs of MTG-Link are the following:
    * Output GFA file: `insertionSequenceCoordinates_mtglink.gfa`
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: `insertionSequenceCoordinates.assembled_sequences.fasta`
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The outputs of MTG-Link are detailed [here](../input-output_files.md)

