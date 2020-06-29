## fasta2gfa.py

The **fasta2gfa.py** script takes as input a FASTA file and treats the 'Ns' regions within the sequences as gaps. It thus transform a FASTA file with 'Ns' regions within sequences to a GFA file.  

It is possible to filter the 'Ns' regions you want to treat as gaps by:
* their size (e.g. gap sizes)
* the flanking contigs' sizes (for example, select only to the 'Ns' regions whose flanking contigs' sizes are long enough to get enough barcodes)

### Usage

```
./fasta2gfa.py --help

usage: fasta2gfa.py -in <fasta_file> -out <output_directory> [options]

Transform a FASTA file with sequences containing 'Ns' regions to a GFA file ('Ns' regions are treated as gaps). We can filter the 'Ns' regions by their size (e.g. gap sizes) and by the contigs' sizes on both sides (long enough for ex to get enough barcodes)
                                
optional arguments:
  -h, --help            show this help message and exit

[Main options]:
  -in INPUT             FASTA file containing the sequences of the scaffolds 
                        obtained from the assembly (format: 'xxx.fasta')
  -min MIN              Minimum size of the 'Ns' region to treat/process as a gap
  -max MAX              Maximum size of the 'Ns' region to treat/process as a gap
  -contigs CONTIGS_SIZE
                        Minimum size of the flanking contigs of the 'Ns' region 
                        to treat/process as a gap
  -out OUTDIR           Output directory for saving the GFA file
```


## paths2gfa.py

The **paths2gfa.py** script takes as input a file containing the paths between scaffolds and transform it to a GFA file.

Format of a path: `<int:nb_scaffolds>****<sid1(f|r)>+<sid2(f|r)>`

### Usage

```
./paths2gfa.py --help

usage: paths2gfa.py -in <fasta_file> -paths <paths_file> -out <output_directory>

Transform a file containing the paths between scaffolds to a GFA file
                                
optional arguments:
  -h, --help            show this help message and exit

[Main options]:
  -in INPUT             FASTA file containing the sequences of the scaffolds 
                        obtained from the assembly (format: 'xxx.fasta')
  -paths PATHS          File containing the paths between scaffolds (obtained from 
                        the matrix) (format: 'xxx.paths.txt')
  -out OUTDIR           Output directory for saving the GFA file and the 
                        corresponding FASTA file
```


## gfa2_to_gfa1.py

The **gfa2_to_gfa1.py** script converts a GFA 2.0 file into a GFA 1.0 file.

### Usage

```
./gfa2_to_gfa1.py --help

usage: gfa.2_to_gfa.1.py -in <input_gfa_2.0)> -out <output_directory>

Convert a GFA 2.0 file into a GFA 1.0 file                                

optional arguments:
  -h, --help            show this help message and exit

[Main options]:
  -in INPUT             GFA 2.0 file (format: 'xxx.gfa')
  -out OUTDIR           Output directory for saving the GFA 1.0 file
```


## gfa2fasta.py

The **gfa2fasta** script takes as input a GFA file (GFA 1.0) and transform it to a FASTA file. The gaps of the GFA file are returned as 'Ns' regions in the output FASTA file.

### Usage

```
./gfa2fasta.py --help

usage: gfa2fasta.py -in <gfa_file> -out <output_directory>

Transform a GFA file (GFA 1.0) to a FASTA file (gaps are returned as 'Ns' regions)
                                
optional arguments:
  -h, --help            show this help message and exit

[Main options]:
  -in INPUT             GFA 1.0 file (format: 'xxx.gfa')
  -out OUTDIR           Output directory for saving the FASTA file
```
