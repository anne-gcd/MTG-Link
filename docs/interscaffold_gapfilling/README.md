# Inter-scaffold gap-filling

Inter-scaffold gap-filling is defined as filling the gaps between scaffolds' sequences (in order to connect scaffolds).  
Below is an example of the procedure used to perform inter-scaffold gap-filling with MTG-Link.  
An example of the outputs you should get is located in the `outputs/` directory.


## How to generate the GFA file

**Input file**: A **matrix** file representing the links between scaffolds in tabular format
* the matrix file contains the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds).  
**NB**: The matrix file can be obtained with `LRez compare` using a file containing regions of interest in format `chromosome:startPosition-endPosition` (regionsFile.lst)
```
LRez compare -b bamFile.bam -r regionsFile.lst -o matrixFile.matrix
```

**Scripts**: To generate a GFA file from a matrix file containing the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds), you have to use the script **`matrix2gfa.py`** (located in utils/).
* The `matrix2gfa.py` script takes as input a file containing the number of common barcodes between all possibles pairs of scaffolds' extremities in tabular format (matrix) and converts it to a GFA file (GFA 2.0)
    * it is possible to specify the minimal number of links (minimal number of common barcodes) two scaffolds must share to consider a real link between them and so to add the corresponding gap to the GFA file: `-threshold` option
```
# Create a GFA file from the matrix file, with a minimal number of links between two scaffolds being of 10
../../utils/matrix2gfa.py -fa interscaffold_gapfilling.fasta -matrix interscaffold_gapfilling.matrix -out . -threshold 10
```
* interscaffold_gapfilling.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* interscaffold_gapfilling.matrix: Matrix file containing the number of common barcodes between all possibles pairs of scaffolds' extremities
* **Outputs**: 
    * GFA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa'
    * FASTA files containing the scaffolds sequences: 'interscaffold_gapfilling_[scaffoldID].scaffold.fasta'


## MTG-Link pipeline

```
# Build the LRez barcode index
LRez index fastq -f interscaffold_gapfilling.fastq.gz -o interscaffold_gapfilling.bci -g

# Run MTG-Link
../../mtglink.py DBG -gfa interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa -bam interscaffold_gapfilling.bam -fastq interscaffold_gapfilling.fastq.gz -index interscaffold_gapfilling.bci -t 4 -k 61 51 41 31 21
```
* interscaffold_gapfilling.fastq.gz: Linked reads file (gzipped), with the barcode sequence in the header (BX:Z tag)
* interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.gfa: GFA file obtained from `matrix2gfa.py`
* interscaffold_gapfilling.bam: input BAM file
* interscaffold_gapfilling.bci: LRez barcode index of the FASTQ file, obtained with `LRez index fastq`
* **Outputs**: the outputs are located in the `MTG-Link_results/` directory
    * Output GFA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa'
        * it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequences.
    * Output FASTA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10.assembled_sequences.fasta'
        * it is a sequence file in FASTA format, that contains the set of assembled target sequences.

**NB**: The output files and directory of MTG-Link are detailed [here](TODO)


## Get the updated FASTA file with assembled sequences

The GFA file output by MTG-Link is an assembly graph that represents the relationships between the input sequences and the target sequences (assemblies or gaps). The input FASTA file can be updated with the assembled sequences using this GFA file. The gaps of the GFA file are returned as 'Ns' regions in the output FASTA file.
```
# Get a BED file containing the coordinates of the 'Ns' regions
../../utils/gfa2tofasta.py -in ./MTG-Link_results/interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa -out .
```
* interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink.gfa: GFA file output by ***MTG-Link***
* **Output**:
    * FASTA file: 'interscaffold_gapfilling_interscaffold_gapfilling_threshold_10_mtglink_assembly.fasta'

