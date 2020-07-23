### Test the installation of MTG-Link

Simulated run with solutions obtained for the higher values of `-k` and `-a`:  
`mtglink.py -gfa test.gfa -c 5000 -bam test.bam -fastq reads.sorted.fastq -index barcoded.shelve -out results_MTGLink`

Results:  
For this run, you should get the following:  
* gap 8-L+_8-R+: 2 solutions for k51.a3  
    * 1 forward of length 2000 bp (Quality AAA)  
    * 1 reverse of length 2000 bp (Quality AAA)