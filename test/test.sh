#! /bin/bash

## Test1: Installation of MTG-Link

echo "TEST1: INSTALLATION"

### MTG-Link DBG
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py DBG -gfa test/test_install/test.gfa -bam test/test_install/test.bam -fastq test/test_install/reads.sorted.fastq.gz -index test/test_install/barcoded.bci -out test/test_install/results_MTGLink_DBG > LOGMTGLinkDBGInstall 2> ERRLOGMTGLinkDBGInstall  
if awk '/\t\* 8_0-3152098-L+:8_3153098-6305195-R+\n\t\t\* fwd\t2000 bp\tAA\n\t\t\* rev\t2000 bp\tAA/' LOGMTGLinkDBGInstall
then
    echo "MTG-Link DBG: Pass"
else
    echo "MTG-Link DBG: Fail"
fi

rm LOGMTGLinkDBGInstall
rm ERRLOGMTGLinkDBGInstall
rm -rf test_install/results_MTGLink_DBG

### MTG-Link IRO
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py IRO -gfa test/test_install/test.gfa -bam test/test_install/test.bam -fastq test/test_install/reads.sorted.fastq.gz -index test/test_install/barcoded.bci -out test/test_install/results_MTGLink_IRO > LOGMTGLinkIROInstall 2> ERRLOGMTGLinkIROInstall  
if awk '/\t\* 8_0-3152098-L+:8_3153098-6305195-R+\n\t\t\* fwd\t2000 bp\tAA/' LOGMTGLinkIROInstall
then
    echo "MTG-Link IRO: Pass"
else
    echo "MTG-Link IRO: Fail"
fi

rm LOGMTGLinkIROInstall
rm ERRLOGMTGLinkIROInstall
rm -rf test_install/results_MTGLink_IRO


## Test2: Run of MTG-Link for gaps INTO scaffolds

echo "TEST2: RUN: GAPS INTO SCAFFOLDS"

### MTG-Link DBG without reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py DBG -gfa test/test_run/gapsIntoScaffolds/test.gfa -bam test/test_run/gapsIntoScaffolds/test.bam -fastq test/test_run/gapsIntoScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsIntoScaffolds/barcoded.bci -k 61 51 41 31 21 -t 3 -out test/test_run/gapsIntoScaffolds/results_MTGLink_DBG > LOGMTGLinkDBGRunInto 2> ERRLOGMTGLinkDBGRunInto  
if (awk '/\t\* 68_0-222132-L+:68_223132-445265-R+\n\t\t\* fwd\t2000 bp\tAA\n\t\t\* rev\t2000 bp\tAA\n\t\* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t\* fwd\t2004 bp\tAA\n\t\t\* rev\t2004 bp\tAA/' LOGMTGLinkDBGRunInto) && (awk '/The gap 27292_0-183358-L+_27292_184358-367715-R+ was not successfully gap-filled/' LOGMTGLinkDBGRunInto)
then
    echo "MTG-Link DBG w/o ref: Pass"
else
    echo "MTG-Link DBG w/o ref: Fail"
fi

rm LOGMTGLinkDBGRunInto
rm ERRLOGMTGLinkDBGRunInto
rm -rf test/test_run/gapsIntoScaffolds/results_MTGLink_DBG

### MTG-Link IRO without reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py IRO -gfa test/test_run/gapsIntoScaffolds/test.gfa -bam test/test_run/gapsIntoScaffolds/test.bam -fastq test/test_run/gapsIntoScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsIntoScaffolds/barcoded.bci -s 20 -t 3 -out test/test_run/gapsIntoScaffolds/results_MTGLink_IRO > LOGMTGLinkIRORunInto 2> ERRLOGMTGLinkIRORunInto  
if (awk '/\t\* 68_0-222132-L+:68_223132-445265-R+\n\t\t\* fwd\t1998 bp\tAB\n\t\* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t\* fwd\t2014 bp\tAA/' LOGMTGLinkIRORunInto) && (awk '/The gap 27292_0-183358-L+_27292_184358-367715-R+ was not successfully gap-filled/' LOGMTGLinkIRORunInto)
then
    echo "MTG-Link IRO w/o ref: Pass"
else
    echo "MTG-Link IRO w/o ref: Fail"
fi

rm LOGMTGLinkIRORunInto
rm ERRLOGMTGLinkIRORunInto
rm -rf test/test_run/gapsIntoScaffolds/results_MTGLink_IRO

### MTG-Link DBG with reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py DBG -gfa test/test_run/gapsIntoScaffolds/test.gfa -bam test/test_run/gapsIntoScaffolds/test.bam -fastq test/test_run/gapsIntoScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsIntoScaffolds/barcoded.bci -k 61 51 41 31 21 -t 3 -refDir test/test_run/gapsIntoScaffolds -out test/test_run/gapsIntoScaffolds/results_MTGLink_DBG_withRef > LOGMTGLinkDBGRunIntoRef 2> ERRLOGMTGLinkDBGRunIntoRef  
if (awk '/\t\* 68_0-222132-L+:68_223132-445265-R+\n\t\t\* fwd\t2000 bp\tA\n\t\t\* rev\t2000 bp\tA\n\t\* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t\* fwd\t2004 bp\tB\n\t\t\* rev\t2004 bp\tB/' LOGMTGLinkDBGRunIntoRef) && (awk '/The gap 27292_0-183358-L+_27292_184358-367715-R+ was not successfully gap-filled/' LOGMTGLinkDBGRunIntoRef)
then
    echo "MTG-Link DBG w/ ref: Pass"
else
    echo "MTG-Link DBG w/ ref: Fail"
fi

rm LOGMTGLinkDBGRunIntoRef
rm ERRLOGMTGLinkDBGRunIntoRef
rm -rf test/test_run/gapsIntoScaffolds/results_MTGLink_DBG_withRef

### MTG-Link IRO with reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py IRO -gfa test/test_run/gapsIntoScaffolds/test.gfa -bam test/test_run/gapsIntoScaffolds/test.bam -fastq test/test_run/gapsIntoScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsIntoScaffolds/barcoded.bci -s 20 -t 3 -refDir test/test_run/gapsIntoScaffolds -out test/test_run/gapsIntoScaffolds/results_MTGLink_IRO_withRef > LOGMTGLinkIRORunIntoRef 2> ERRLOGMTGLinkIRORunIntoRef  
if (awk '/\t\* 68_0-222132-L+:68_223132-445265-R+\n\t\t\* fwd\t1998 bp\tB\n\t\* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t\* fwd\t2014 bp\tB/' LOGMTGLinkIRORunIntoRef) && (awk '/The gap 27292_0-183358-L+_27292_184358-367715-R+ was not successfully gap-filled/' LOGMTGLinkIRORunIntoRef)
then
    echo "MTG-Link IRO w/ ref: Pass"
else
    echo "MTG-Link IRO w/ ref: Fail"
fi

rm LOGMTGLinkIRORunIntoRef
rm ERRLOGMTGLinkIRORunIntoRef
rm -rf test/test_run/gapsIntoScaffolds/results_MTGLink_IRO_withRef


## Test3: Run of MTG-Link for gaps BETWEEN scaffolds

echo "TEST3: RUN: GAPS BETWEEN SCAFFOLDS"

### MTG-Link DBG without reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py DBG -gfa test/test_run/gapsBetweenScaffolds/test.gfa -bam test/test_run/gapsBetweenScaffolds/test.bam -fastq test/test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsBetweenScaffolds/barcoded.bci -k 61 51 41 31 21 -t 3 -out test/test_run/gapsBetweenScaffolds/results_MTGLink_DBG > LOGMTGLinkDBGRunBetween 2> ERRLOGMTGLinkDBGRunBetween
if (awk '/\t\* 68-1+:68-2+\n\t\t\* fwd\t2000 bp\tAA\n\t\t\* rev\t2000 bp\tAA\n\t\* 26939-1+:26939-2+\n\t\t\* fwd\t2004 bp\tAA\n\t\t\* rev\t2004 bp\tAA/' LOGMTGLinkDBGRunBetween) && (awk '/The gap 27292-1+_27292-2+ was not successfully gap-filled/' LOGMTGLinkDBGRunBetween)
then
    echo "MTG-Link DBG w/o ref: Pass"
else
    echo "MTG-Link DBG w/o ref: Fail"
fi

rm LOGMTGLinkDBGRunBetween
rm ERRLOGMTGLinkDBGRunBetween
rm -rf test/test_run/gapsBetweenScaffolds/results_MTGLink_DBG

### MTG-Link IRO without reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py IRO -gfa test/test_run/gapsBetweenScaffolds/test.gfa -bam test/test_run/gapsBetweenScaffolds/test.bam -fastq test/test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsBetweenScaffolds/barcoded.bci -s 20 -t 3 -out test/test_run/gapsBetweenScaffolds/results_MTGLink_IRO > LOGMTGLinkIRORunBetween 2> ERRLOGMTGLinkIRORunBetween
if (awk '/\t\* 68-1+:68-2+\n\t\t\* fwd\t1998 bp\tAB\n\t\* 26939-1+:26939-2+\n\t\t\* fwd\t2014 bp\tAA/' LOGMTGLinkIRORunBetween) && (awk '/The gap 27292-1+_27292-2+ was not successfully gap-filled/' LOGMTGLinkIRORunBetween)
then
    echo "MTG-Link IRO w/o ref: Pass"
else
    echo "MTG-Link IRO w/o ref: Fail"
fi

rm LOGMTGLinkIRORunBetween
rm ERRLOGMTGLinkIRORunBetween
rm -rf test/test_run/gapsBetweenScaffolds/results_MTGLink_IRO

### MTG-Link DBG with reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py DBG -gfa test/test_run/gapsBetweenScaffolds/test.gfa -bam test/test_run/gapsBetweenScaffolds/test.bam -fastq test/test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsBetweenScaffolds/barcoded.bci -k 61 51 41 31 21 -t 3 -refDir test/test_run/gapsBetweenScaffolds -out test/test_run/gapsBetweenScaffolds/results_MTGLink_DBG_withRef > LOGMTGLinkDBGRunBetweenRef 2> ERRLOGMTGLinkDBGRunBetweenRef  
if (awk '/\t\* 68-1+:68-2+\n\t\t\* fwd\t2000 bp\tA\n\t\t\* rev\t2000 bp\tA\n\t\* 26939-1+:26939-2+\n\t\t\* fwd\t2004 bp\tB\n\t\t\* rev\t2004 bp\tB/' LOGMTGLinkDBGRunBetweenRef) && (awk '/The gap 27292-1+_27292-2+ was not successfully gap-filled/' LOGMTGLinkDBGRunBetweenRef)
then
    echo "MTG-Link DBG w/ ref: Pass"
else
    echo "MTG-Link DBG w/ ref: Fail"
fi

rm LOGMTGLinkDBGRunBetweenRef
rm ERRLOGMTGLinkDBGRunBetweenRef
rm -rf test/test_run/gapsBetweenScaffolds/results_MTGLink_DBG_withRef

### MTG-Link IRO with reference
. /local/env/envconda.sh
conda activate mtglink

/usr/bin/time -v ../mtglink.py IRO -gfa test/test_run/gapsBetweenScaffolds/test.gfa -bam test/test_run/gapsBetweenScaffolds/test.bam -fastq test/test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz -index test/test_run/gapsBetweenScaffolds/barcoded.bci -s 20 -t 3 -refDir test/test_run/gapsBetweenScaffolds -out test/test_run/gapsBetweenScaffolds/results_MTGLink_IRO_withRef > LOGMTGLinkIRORunBetweenRef 2> ERRLOGMTGLinkIRORunBetweenRef  
if (awk '/\t\* 68-1+:68-2+\n\t\t\* fwd\t1998 bp\tB\n\t\* 26939-1+:26939-2+\n\t\t\* fwd\t2014 bp\tB/' LOGMTGLinkIRORunBetweenRef) && (awk '/The gap 27292-1+_27292-2+ was not successfully gap-filled/' LOGMTGLinkIRORunBetweenRef)
then
    echo "MTG-Link IRO w/ ref: Pass"
else
    echo "MTG-Link IRO w/ ref: Fail"
fi

rm LOGMTGLinkIRORunBetweenRef
rm ERRLOGMTGLinkIRORunBetweenRef
rm -rf test/test_run/gapsBetweenScaffolds/results_MTGLink_IRO_withRef

