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
    * the parameter `-extension` allows to define a gap larger than the single position of the insertion event (recommended if the insertion site position is uncertain or not exactly defined). The default value is 50 bp at each side of the position given in the VCF.
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
* the parameter `-m` is set here to 400 bp in order to report only assembled sequences that are larger than 400 bp, that is corresponding to insertion variants that are larger than 100 bp. The assembled sequence contains 2 flanking sequences of 100 bp (`-ext` value, for qualitative evaluation), the insertion gap of 100 bp (defined in the gfa file if it was obtained with `vcf2gfa.py` and its `-extension` parameter set to 50 bp at each side of the insertion site position) and the insertion variant sequence itself. Thus, minimal insertion size = `-m` -2* `-ext` -2* `-extension`
* **Outputs**: the main outputs of MTG-Link are the following:
    * Output GFA file: `insertionSequenceCoordinates_mtglink.gfa`
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: `insertionSequenceCoordinates.assembled_sequences.fasta`
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The outputs of MTG-Link are detailed [here](../input-output_files.md). 

**NB2**: The outputs of MTG-Link can be converted to VCF format using the external script [create_mtglink_vcf.py](https://gist.github.com/pontushojer/da2c40de4f8d89c23fa992d6ff6f3cc3) (author: [Pontus Höjer](https://github.com/pontushojer), see also the discussion in [issue #25](https://github.com/anne-gcd/MTG-Link/issues/25))

