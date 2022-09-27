# Inter-scaffold gap-filling

Inter-scaffold gap-filling is defined as filling the gaps between scaffolds' sequences (in order to connect scaffolds).  
Below is an example of the procedure used to perform inter-scaffold gap-filling with MTG-Link.  
In this example, we are interested to fill gaps between different scaffolds. In this use case, we do not know the order and orientation of the scaffolds of interest, but we use shared barcode links between them to identify putative pairs of consecutive scaffolds (defining the gaps to assemble).  
An example of the outputs you should get is located in the `outputs/` directory.


## How to generate the GFA file

**Input file**: 
* A **matrix** file representing the links between scaffolds in tabular format
    * the matrix file contains the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds).  
    **NB**: The matrix file can be obtained with `LRez compare` using a set of scaffolds of interest (more precisely, a file containing regions of interest in format `chromosome:startPosition-endPosition` (regionsFile.lst)
```
LRez compare -b bamFile.bam -r regionsFile.lst -o matrixFile.matrix
```

**Scripts**: To generate a GFA file from a matrix file containing the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds), you have to use the script **`matrix2gfa.py`** (located in `utils/`).
* The `matrix2gfa.py` script takes as input a file containing the number of common barcodes between all possibles pairs of scaffolds' extremities in tabular format (matrix) and converts it to a GFA file (GFA 2.0)
    * it is possible to specify the minimal number of links (minimal number of common barcodes) two scaffolds must share to consider a real link between them and so to add the corresponding gap to the GFA file: `-threshold` option
```
# Create a GFA file from the matrix file, with a minimal number of links between two scaffolds being of 10
../../utils/matrix2gfa.py -fa referenceGenome.fasta -matrix matrixFile.matrix -out . -threshold 10
```
* referenceGenome.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* matrixFile.matrix: Matrix file containing the number of common barcodes between all possibles pairs of scaffolds' extremities
* **Outputs**: 
    * GFA file: 'referenceGenome_matrixFile_threshold_10.gfa'
    * FASTA files containing the scaffolds sequences: 'referenceGenome_[scaffoldID].scaffold.fasta'


## MTG-Link pipeline

```
# Build the LRez barcode index (if not already done)
LRez index fastq -f readsFile.fastq.gz -o barcodeIndex.bci -g

# Run MTG-Link
../../mtglink.py DBG -gfa referenceGenome_matrixFile_threshold_10.gfa -bam bamFile.bam -fastq readsFile.fastq.gz -index barcodeIndex.bci -t 4 -k 61 51 41 31 21
```
* readsFile.fastq.gz: Linked-reads file (must be gzipped), with the barcode sequence in the header (BX:Z tag)
* referenceGenome_matrixFile_threshold_10.gfa: GFA file obtained from `matrix2gfa.py`
* bamFile.bam: input BAM file
* barcodeIndex.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the main outputs of MTG-Link are the following:
    * Output GFA file: 'referenceGenome_matrixFile_threshold_10_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'referenceGenome_matrixFile_threshold_10.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The outputs of MTG-Link are detailed [here](../input-output_files.md)


## Get the updated FASTA file with assembled sequences

The GFA file output by MTG-Link is an assembly graph that represents the relationships between the input sequences and the target sequences (assemblies or gaps). The input FASTA file can be updated with the assembled sequences using this GFA file. The gaps of the GFA file are returned as 'N's regions in the output FASTA file.
```
# Update the input FASTA file with assembled sequences
../../utils/gfa2tofasta.py -in ./MTG-Link_results/referenceGenome_matrixFile_threshold_10_mtglink.gfa -out .
```
* interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa: GFA file output by ***MTG-Link***
* **Output**:
    * FASTA file: 'referenceGenome_matrixFile_threshold_10_mtglink_assembly.fasta'

