## Test the installation of MTG-Link

### Module DBG

Simulated run with default values:  
`mtglink.py DBG -gfa test.gfa -bam test.bam -fastq reads.sorted.fastq.gz -index barcoded.shelve -out results_MTGLink_DBG`

Results:  
For this run, you should get the following:  
* gap '8_gap1-L+:8_gap1-R+': 2 solutions
    * 1 forward of length 2000 bp (Quality A)  
    * 1 reverse of length 2000 bp (Quality A)

### Module IRO

Simulated run with default values:  
`mtglink.py IRO -gfa test.gfa -bam test.bam -fastq reads.sorted.fastq.gz -index barcoded.shelve -out results_MTGLink_IRO`

Results:  
For this run, you should get the following:  
* gap '8_gap1-L+:8_gap1-R+': 1 solution
    * 1 forward of length 2000 bp (Quality A)