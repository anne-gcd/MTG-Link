# MTG10X

## About MTG10X

The **MTG10X** tool uses other tools relying on GATB-CORE library.

To run the mtg10x.py script, you need to install a virtual environment containing Samtools, Biopython, Gfapy, Blast and Mummer librairies.  
You also need to install the Con10X and MindTheGap tools and add their locations to the PATH:

    * Con10X: <https://github.com/flegeai/Con10x>
    * MindTheGap: <https://github.com/GATB/MindTheGap>

We are using the tools BamExtractor and GetReads, from the **Con10x** repository.  
**MindTheGap** is used in *'breakpoint'* mode, by removing first a small region on both sides (e.g. extension of the gap) of size *-ext* (determines start/end of gapfilling).  
MindTheGap will try to find a path in the **de Bruijn graph** from the left k-mer (source) to the right k-mer (target), and it will try for several k-mer sizes and several abundance thresholds for solid k-mers. However, once it has find a path (e.g. a solution), it will stop searching for other k-mer sizes and abundance thresholds.

When a reference sequence is provided for the gap, some statistics are performed on the results.

Command-line to execute **mtg10x.py**:  
`./mtg10x.py -in <GFA_file> -c <chunk_size> -bam <BAM_file> -reads <reads_file> -index <index_file> [options]
`

Parameters:
    * '-c': chunk size
    * '-f': minimal frequence of barcodes in union (default 2)
    * '-k': size of k-mers (default [51, 41, 31, 21])
    * '-a': minimal abundance threshold for solid k-mers (default [3, 2])
    * '-ext': size of the extension of the gap, on both sides (default k); determines start/end of gapfilling
    * '--force': force search on all k-mer sizes
    * '-out': output directory (default './mtg10x_results')
    * '-refDir': directory where the reference sequence(s) is/are (if any)


## License

Please not that GATB-Core is distributed under Affero-GPL license.

## Dependencies

The following third parties should be already installed:
```
    cmake 3.1+ *[mandatory]*
```

## Packages

The list of packages to install are in requirements.txt
To install a list of packages into a specified conda environment, do the following:
```
    conda create --name <env> --file requirements.txt
```