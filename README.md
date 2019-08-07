#About MTG10X

The MTG10X tool uses other tools relying on GATB-CORE library.

To run the mtg10x.py script, you need to install a virtual environment containing Samtools, Biopython and Gfapy librairies.
You also need to pre-install the Con10X and MindTheGap tools on your working directory:
'''
    * Con10X: https://github.com/flegeai/Con10x
    * MindTheGap: https://github.com/GATB/MindTheGap
'''
We are using the tools BamExtractor and GetReads, from the Con10x repository.
MindTheGap is used in 'breakpoint' mode, and with removed flanking regions of size k.

This tool will be further completed with statistics on the results (ongoing).

Command-line to execute mtg10x.py:
'''
    ./mtg10x.py -gfa <GFA_file> -c <chunk_size> -bam <BAM_file> -reads <reads_file> -index <index_file> [options]
'''

#License

Please not that GATB-Core is distributed under Affero-GPL license.

#Dependencies

The following third parties should be already installed:
'''
    * cmake 3.1+ (mandatory)
'''

#Packages

The list of packages to install are in requirements.txt
To install a list of packages into a specified conda environment, do the following:
'''
    conda create --name <env> --file requirements.txt
'''