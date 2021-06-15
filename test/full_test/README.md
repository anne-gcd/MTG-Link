### Tests on a simulated run

#### Module DBG

To test MTG-Link for multiple features:  
* managing the orientations '+' and '-' of the flanking scaffolds of the gap  
* automatic gap-filling with decreasing values of `-k` and `-a`  
* qualitative evaluation performed with the flanking contigs information  
* qualitative evaluation performed with reference sequences

##### Simulated dataset without reference sequences

Simulated run with solutions obtained for different values of `-k` and `-a`, the qualitative evaluation being performed with the flanking contigs information:  
`mtglink.py DBG -gfa test.gfa -bam test.bam -fastq reads.sorted.fastq -index barcoded.shelve -k 61 51 41 31 21 -out results_MTGLink`  

Results:  
For this run, you should get the following:  
* gap '95-L-_95-R-': 2 solutions
    * 1 forward of length 2000 bp (Quality AA)  
    * 1 reverse of length 2000 bp (Quality AA)  
* gap '68-L+_68-R+': 2 solutions 
    * 1 forward of length 2000 bp (Quality AA)  
    * 1 reverse of length 2000 bp (Quality AA)  
* gap '26939-L+_26939-R+': 2 solutions
    * 1 forward of length 2004 bp (Quality AA)  
    * 1 reverse of length 2004 bp (Quality AA)  
* gap '27292-L+_27292-R+': no solution

##### Simulated dataset with reference sequences

Simulated run with solutions obtained for different values of `-k` and `-a`, the qualitative evaluation being performed with the reference sequence:  
`mtglink.py DBG -gfa test.gfa -bam test.bam -fastq reads.sorted.fastq -index barcoded.shelve -k 61 51 41 31 21 -refDir . -out results_MTGLink_withref`

Results:  
For this run, you should get the following:  
* gap '95-L-_95-R-': 2 solutions
    * 1 forward of length 2000 bp (Quality A)  
    * 1 reverse of length 2000 bp (Quality A)  
* gap '68-L+_68-R+': 2 solutions  
    * 1 forward of length 2000 bp (Quality A)  
    * 1 reverse of length 2000 bp (Quality A)  
* gap '26939-L+_26939-R+': 2 solutions
    * 1 reverse of length 2004 bp (Quality B)  
    * 1 forward of length 2004 bp (Quality B)  
* gap '27292-L+_27292-R+': no solution

