# Get a GFA file from a FASTA file containing 'Ns' regions

## fasta2bed.py

The **`fasta2bed.py`** script takes as input a FASTA file and converts it to a BED file containing the 'N's coordinates for each scaffold. It can be run with the following command:  
```
fasta2bed.py -fa referenceGenome.fasta -out targetSequenceCoordinates.bed
```
* referenceGenome.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* targetSequenceCoordinates.bed: Output BED file

The format of the output BED file is: `<scaffoldID>  <Nstretch_startPosition>  <Nstretch_endPosition>`

## bed2gfa.py

The **`bed2gfa.py`** script converts a BED file containing the 'N's coordinates for each scaffold (or locus coordinates) to a GFA file (GFA 2.0).  
It is possible to filter the 'N's regions you want to treat as gaps by:
* their size (e.g. gap lengths)
* the flanking contigs' sizes (for example, select only the 'N's regions whose flanking contigs' sizes are long enough to get enough barcodes)

It can be run with the following command:  
```
bed2gfa.py -bed targetSequenceCoordinates.bed -fa referenceGenome.fasta -out targetSequenceCoordinates.gfa -min MIN_GAPLENGTH -max MAX_GAPLENGTH -contigs MIN_CONTIGSIZE
```
* targetSequenceCoordinates.bed: BED file containing the 'Ns' coordinates for each scaffold
* referenceGenome.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* targetSequenceCoordinates.gfa: Output GFA file
* MIN_GAPLENGTH: Minimum size of the 'N's region to treat as a gap
* MAX_GAPLENGTH: Maximum size of the 'N's region to treat as a gap
* MIN_CONTIGSIZE: Minimum size of the flanking contigs of the 'N's region to treat as a gap

The output GFA file is a GFA 2.0 file. 


# Get a GFA file from a file containing the paths between scaffolds

## paths2gfa.py

The **`paths2gfa.py`** script takes as input a file containing the paths between scaffolds and converts it to a GFA file (GFA 2.0). It can be run with the following command:  
```
paths2gfa.py -fa referenceGenome.fasta -paths pathsFile.txt -out targetSequenceCoordinates.gfa
```
* referenceGenome.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* pathsFile.txt: File containing the paths between scaffolds
* targetSequenceCoordinates.gfa: Output GFA file

The format of the input PATHS file is: `<numberOfScaffolds>****<sid1(f|r)>+<sid2(f|r)>`  
When the orientation of the scaffold is undetermined (`?`), both forward and reverse orientations are taken into consideration.  
The output GFA file is a GFA 2.0 file. 


# Get a GFA file from a matrix file containing the links between the ends of the scaffolds

## matrix2gfa.py

The **`matrix2gfa.py`** script takes as input a file containing the number of common barcodes between all possibles pairs of scaffolds' extremities (links between the ends of the scaffolds) in tabular format (matrix) and converts it to a GFA file (GFA 2.0). It can be run with the following command:  
```
matrix2gfa.py -fa referenceGenome.fasta -matrix matrixFile.matrix -out targetSequenceCoordinates.gfa -threshold THRESHOLD
```
* referenceGenome.fasta: FASTA file containing the sequences of the scaffolds (reference genome)
* matrixFile.matrix: File containing the number of common barcodes between all possibles pairs of scaffolds' extremities (matrix)
* targetSequenceCoordinates.gfa: Output GFA file
* THRESHOLD: Minimal number of links two scaffolds must share to try to fill the gap between them

The format of the input MATRIX file is: `<scaffoldLeftID>:<extremitiesCoordStart>-<extremitiesCoordEnd> <scaffoldRightID>:<extremitiesCoordStart>-<extremitiesCoordEnd> <numberOfCommonBarcodes>`  
The output GFA file is a GFA 2.0 file.


# Get a GFA file from a VCF file containing Indels coordinates

## vcf2gfa.py

The **`vcf2gfa.py`** script converts a VCF file containing the Indels coordinates to a GFA file (GFA 2.0).  
It is possible to filter the Indels you want to treat as targets (e.g. you want to reconstruct) by:
* the flanking contigs' sizes (for example, select only the Indels whose flanking contigs' sizes are long enough to get enough barcodes)

It can be run with the following command:  
```
vcf2gfa.py -vcf vcfFile.vcf -fa referenceGenome.fasta -out indelSequenceCoordinates.gfa -contigs MIN_CONTIGSIZE
```
* vcfFile.vcf: VCF file containing the Indels coordinates
* referenceGenome.fasta: FASTA file of the reference genome
* indelSequenceCoordinates.gfa: Output GFA file
* MIN_CONTIGSIZE: Minimum size of the flanking contigs of the Indels to treat as a target

The output GFA file is a GFA 2.0 file. 


# Merge two GFA 2.0 files together

## mergegfa.py

The **mergegfa.py** script takes as input two GFA files (GFA 2.0) and merge them together. It can be run with the following command:  
```
mergegfa.py -1 gfaFile1.gfa -2 gfaFile2.gfa -out merged_gfa.gfa
```
* gfaFile1.gfa: GFA 2.0 file n°1
* gfaFile2.gfa: GFA 2.0 file n°2
* merged_gfa.gfa: Output merged GFA file

The output GFA file is a GFA 2.0 file. 


# Convert a GFA 2.0 file to a FASTA file

## gfa2tofasta.py

The **`gfa2tofasta.py`** script takes as input a GFA file (GFA 2.0) and converts it to a FASTA file. The gaps of the GFA file are returned as 'N's regions in the output FASTA file. It can be run with the following command:  
```
gfa2tofasta.py -in gfaFile.gfa -out fastaFile.fasta
```
* gfaFile.gfa: GFA 2.0 file
* fastaFile.fasta: Output FASTA file

