### Tests on a simulated run

To test MTG-Link for multiple features:  
* managing the orientations '+' and '-' of the flanking scaffolds of the gap  
* automatic gap-filling with decreasing values of `-k` and `-a`  
* qualitative evaluation performed with the flanking contigs information  
* qualitative evaluation performed with reference sequences

#### Simulated dataset without reference sequences

Simulated run with solutions obtained for different values of `-k` and `-a`, the qualitative evaluation being performed with the flanking contigs information:  
`mtglink.py -gfa test.gfa -c 5000 -bam test.bam -fastq reads.sorted.fastq -index barcoded.shelve -k 61 51 41 31 21 -out results_MTGLink`  

Results:  
For this run, you should get the following:  
* gap 95-L-_95-R-: 2 solutions for k61.a2  
    * 1 forward of length 2000 bp (Quality AAA)  
    * 1 reverse of length 2000 bp (Quality AAA)  
* gap 68-L+_68-R+: 2 solutions for k61.a2  
    * 1 forward of length 2000 bp (Quality AAA)  
    * 1 reverse of length 2000 bp (Quality AAA)  
* gap 26939-L+_26939-R+: 2 solution for k51.a3  
    * 1 forward of length 2004 bp (Quality AAA)  
    * 1 reverse of length 2004 bp (Quality AAA)  
* scaffold 27292: no solution

#### Simulated dataset with reference sequences

Simulated run with solutions obtained for different values of `-k` and `-a`, the qualitative evaluation being performed with the reference sequence:  
`mtglink.py -gfa test.gfa -c 5000 -bam test.bam -fastq reads.sorted.fastq -index barcoded.shelve -k 61 51 41 31 21 -refDir test/full_test/ -out results_MTGLink_withref`

Results:  
For this run, you should get the following:  
* gap 95-L-_95-R-: 2 solutions for k61.a2  
    * 1 forward of length 2000 bp (Quality AA)  
    * 1 reverse of length 2000 bp (Quality AA)  
* gap 68-L+_68-R+: 2 solutions for k61.a2  
    * 1 forward of length 2000 bp (Quality AA)  
    * 1 reverse of length 2000 bp (Quality AA)  
* gap 26939-L+_26939-R+: 2 solution for k51.a3  
    * 1 forward of length 2004 bp (Quality BA)  
    * 1 reverse of length 2004 bp (Quality BA)  
* scaffold 27292: no solution

