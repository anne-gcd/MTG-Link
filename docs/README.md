# Examples of MTG-Link applications

MTG-Link can be used for various local assembly use cases, such as the reconstruction of loci of interest, intra-scaffold and inter-scaffold gap-fillings, as well as alternative allele reconstruction of large insertion variants.  

We provide here tiny toy datasets for each of these use cases with all the command lines to run. The main differences between these use cases lie in how to generate the input GFA file, and also how to convert the output of MTG-Link in other formats.

 
## A single locus of interest

We are interested here in a particular locus on a reference genome. We have re-sequenced an individual from the same species with a linked-read technology and we want to assemble this specific locus from the sequencing reads of this individual. 

The locus of interest is defined by two coordinates on a chromosome or scaffold of the reference genome, therefore we dispose of a BED file (locus coordinates) and a sequence file (reference genome) to generate the input GFA file for MTG-Link.

```
# Generate the input GFA file
../utils/bed2gfa.py -bed myLocus.bed -fa referenceGenome.fasta -out myLocus.gfa

# Build the LRez barcode index (if not already done)
LRez index fastq -f readsFile.fastq.gz -o barcodeIndex.bci -g

# Run MTG-Link
../mtglink.py DBG -gfa myLocus.gfa -bam bamFile.bam -fastq readsFile.fastq.gz -index barcodeIndex.bci -k 61 51 41 31 21
```

**Outputs**: the main outputs of MTG-Link are the following:
* Output GFA file: `myLocus_mtglink.gfa`
	* it is an assembly graph file in GFA format, that complements the input GFA file with the obtained assembled target sequence.
* Output FASTA file: `myLocus.assembled_sequences.fasta`
	* it is a sequence file in FASTA format, that contains the assembled target sequence.

**NB**: The outputs of MTG-Link are detailed [here](../input-output_files.md)  
**Note**: if we have not one but several loci to assemble, they can be all given in the same input BED file with one locus per line.

The example of this use case with a tiny dataset is detailed [here](./singlelocus_reconstruction/README.md).


## Intra-scaffold gap-filling

In intra-scaffold gap-filling, the target sequences to be assembled are defined by stretches of 'N's in a sequence file (reference genome). We provide scripts to select the gaps to assemble and to generate the corresponding input GFA file. In addition, for this use case, the input sequence file can be updated at the end of MTG-Link run by replacing the Ns by the assembled sequences.

The full example of this use case with a tiny dataset is detailed [here](./intrascaffold_gapfilling/README.md).


## Inter-scaffold gap-filling

In inter-scaffold gap-filling, we are interested to fill gaps between different scaffolds. In this use case, we do not know the order and orientation of the scaffolds of interest, but we use shared barcode links between them to identify putative pairs of consecutive scaffolds (defining the gaps to assemble). We provide scripts to generate the corresponding input GFA file, starting from a set of scaffolds of interest. In addition, for this use case, the input scaffold file can be updated at the end of MTG-Link run by replacing scaffolded and gap-filled scaffolds with super-scaffolds.

The full example of this use case with a tiny dataset is detailed [here](./interscaffold_gapfilling/README.md).


## Alternative allele reconstruction of large insertion variants

In this application, we have identified insertion variants from a re-sequenced individual with respect to a reference genome, but we dispose only of the predicted insertion site and the inserted sequence is unknown. Insertion sites are defined in a VCF file. We provide scripts to generate the corresponding input GFA file, starting from the VCF file and the reference genome. 

The full example of this use case with a tiny dataset is detailed [here](./insertionsvariants_reconstruction/README.md).