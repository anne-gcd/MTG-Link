#!/usr/bin/env python3

import subprocess


## Test1: Installation of MTG-Link

print("TEST1: INSTALLATION")

### MTG-Link DBG
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_install/test.gfa", "-bam", "test_install/test.bam", \
            "-fastq", "test_install/reads.sorted.fastq.gz", "-index", "test_install/barcoded.bci", \
            "-out", "test_install/results_MTGLink_DBG"]
commandOut = "LOGMTGLinkDBGInstall"
commandLog = "ERRLOGMTGLinkDBGInstall"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 1 targets were successfully assembled:\n\t* 8_0-3152098-L+:8_3153098-6305195-R+\n\t\t* fwd_1/1.k51\t2000 bp\tAA\n\t\t* rev_1/1.k51\t2000 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link DBG: Pass")
    else:
        print("MTG-Link DBG: Fail")

subprocess.run("rm LOGMTGLinkDBGInstall", shell=True)
subprocess.run("rm ERRLOGMTGLinkDBGInstall", shell=True)
subprocess.run("rm -rf test_install/results_MTGLink_DBG", shell=True)


### MTG-Link IRO
command = ["/usr/bin/time", "-v", "../mtglink.py", "IRO", \
            "-gfa", "test_install/test.gfa", "-bam", "test_install/test.bam", \
            "-fastq", "test_install/reads.sorted.fastq.gz", "-index", "test_install/barcoded.bci", \
            "-out", "test_install/results_MTGLink_IRO"]
commandOut = "LOGMTGLinkIROInstall"
commandLog = "ERRLOGMTGLinkIROInstall"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 1 targets were successfully assembled:\n\t* 8_0-3152098-L+:8_3153098-6305195-R+\n\t\t* fwd.k0\t2000 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link IRO: Pass")
    else:
        print("MTG-Link IRO: Fail")

subprocess.run("rm LOGMTGLinkIROInstall", shell=True)
subprocess.run("rm ERRLOGMTGLinkIROInstall", shell=True)
subprocess.run("rm -rf test_install/results_MTGLink_IRO", shell=True)



## Test2: Run of MTG-Link for gaps INTO scaffolds

print("TEST2: RUN: GAPS INTO SCAFFOLDS")

### MTG-Link DBG without reference
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "31", "21", "-t", "3", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG"]
commandOut = "LOGMTGLinkDBGRunInto"
commandLog = "ERRLOGMTGLinkDBGRunInto"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link DBG w/o ref: Pass")
    else:
        print("MTG-Link DBG w/o ref: Fail")

subprocess.run("rm LOGMTGLinkDBGRunInto", shell=True)
subprocess.run("rm ERRLOGMTGLinkDBGRunInto", shell=True)
subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG", shell=True)


### MTG-Link IRO without reference
command = ["/usr/bin/time", "-v", "../mtglink.py", "IRO", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-s", "20", "-t", "3", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_IRO"]
commandOut = "LOGMTGLinkIRORunInto"
commandLog = "ERRLOGMTGLinkIRORunInto"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd.k0\t1998 bp\tAB\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd.k0\t2014 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link IRO w/o ref: Pass")
    else:
        print("MTG-Link IRO w/o ref: Fail")

subprocess.run("rm LOGMTGLinkIRORunInto", shell=True)
subprocess.run("rm ERRLOGMTGLinkIRORunInto", shell=True)
subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_IRO", shell=True)


# ### MTG-Link DBG with reference
# /usr/bin/time -v ../mtglink.py DBG -gfa test_run/gapsIntoScaffolds/test.gfa -bam test_run/gapsIntoScaffolds/test.bam -fastq test_run/gapsIntoScaffolds/reads.sorted.fastq.gz -index test_run/gapsIntoScaffolds/barcoded.bci -k 61 51 41 31 21 -t 3 -refDir test_run/gapsIntoScaffolds -out test_run/gapsIntoScaffolds/results_MTGLink_DBG_withRef > LOGMTGLinkDBGRunIntoRef 2> ERRLOGMTGLinkDBGRunIntoRef  
# val=$(awk '/\t\* 68_0-222132-L+:68_223132-445265-R+\n\t\t\* fwd\t2000 bp\tA\n\t\t\* rev\t2000 bp\tA\n\t\* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t\* fwd\t2004 bp\tB\n\t\t\* rev\t2004 bp\tB/' LOGMTGLinkDBGRunIntoRef)
# val1=$(awk '/The gap 27292_0-183358-L+_27292_184358-367715-R+ was not successfully gap-filled/' LOGMTGLinkDBGRunIntoRef)
# if [ -z "$val" ] || [ -z "$val1" ]
# then
#     echo "MTG-Link DBG w/ ref: Fail"
# else
#     echo "MTG-Link DBG w/ ref: Pass"
# fi

# rm LOGMTGLinkDBGRunIntoRef
# rm ERRLOGMTGLinkDBGRunIntoRef
# rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_withRef

# ### MTG-Link IRO with reference
# /usr/bin/time -v ../mtglink.py IRO -gfa test_run/gapsIntoScaffolds/test.gfa -bam test_run/gapsIntoScaffolds/test.bam -fastq test_run/gapsIntoScaffolds/reads.sorted.fastq.gz -index test_run/gapsIntoScaffolds/barcoded.bci -s 20 -t 3 -refDir test_run/gapsIntoScaffolds -out test_run/gapsIntoScaffolds/results_MTGLink_IRO_withRef > LOGMTGLinkIRORunIntoRef 2> ERRLOGMTGLinkIRORunIntoRef  
# if (awk '/\t\* 68_0-222132-L+:68_223132-445265-R+\n\t\t\* fwd\t1998 bp\tB\n\t\* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t\* fwd\t2014 bp\tB/' LOGMTGLinkIRORunIntoRef) && (awk '/The gap 27292_0-183358-L+_27292_184358-367715-R+ was not successfully gap-filled/' LOGMTGLinkIRORunIntoRef)
# then
#     echo "MTG-Link IRO w/ ref: Pass"
# else
#     echo "MTG-Link IRO w/ ref: Fail"
# fi

# rm LOGMTGLinkIRORunIntoRef
# rm ERRLOGMTGLinkIRORunIntoRef
# rm -rf test_run/gapsIntoScaffolds/results_MTGLink_IRO_withRef



## Test3: Run of MTG-Link for gaps BETWEEN scaffolds

print("TEST3: RUN: GAPS BETWEEN SCAFFOLDS")

### MTG-Link DBG without reference
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsBetweenScaffolds/test.gfa", "-bam", "test_run/gapsBetweenScaffolds/test.bam", \
            "-fastq", "test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsBetweenScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "31", "21", "-t", "3", \
            "-out", "test_run/gapsBetweenScaffolds/results_MTGLink_DBG"]
commandOut = "LOGMTGLinkDBGRunBetween"
commandLog = "ERRLOGMTGLinkDBGRunBetween"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68-1+:68-2+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t* 26939-1+:26939-2+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link DBG w/o ref: Pass")
    else:
        print("MTG-Link DBG w/o ref: Fail")

subprocess.run("rm LOGMTGLinkDBGRunBetween", shell=True)
subprocess.run("rm ERRLOGMTGLinkDBGRunBetween", shell=True)
subprocess.run("rm -rf test_run/gapsBetweenScaffolds/results_MTGLink_DBG", shell=True)


### MTG-Link IRO without reference
command = ["/usr/bin/time", "-v", "../mtglink.py", "IRO", \
            "-gfa", "test_run/gapsBetweenScaffolds/test.gfa", "-bam", "test_run/gapsBetweenScaffolds/test.bam", \
            "-fastq", "test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsBetweenScaffolds/barcoded.bci", \
            "-s", "20", "-t", "3", \
            "-out", "test_run/gapsBetweenScaffolds/results_MTGLink_IRO"]
commandOut = "LOGMTGLinkIRORunBetween"
commandLog = "ERRLOGMTGLinkIRORunBetween"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68-1+:68-2+\n\t\t* fwd.k0\t1998 bp\tAB\n\t* 26939-1+:26939-2+\n\t\t* fwd.k0\t2014 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link IRO w/o ref: Pass")
    else:
        print("MTG-Link IRO w/o ref: Fail")

subprocess.run("rm LOGMTGLinkIRORunBetween", shell=True)
subprocess.run("rm ERRLOGMTGLinkIRORunBetween", shell=True)
subprocess.run("rm -rf test_run/gapsBetweenScaffolds/results_MTGLink_IRO", shell=True)


# ### MTG-Link DBG with reference
# /usr/bin/time -v ../mtglink.py DBG -gfa test_run/gapsBetweenScaffolds/test.gfa -bam test_run/gapsBetweenScaffolds/test.bam -fastq test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz -index test_run/gapsBetweenScaffolds/barcoded.bci -k 61 51 41 31 21 -t 3 -refDir test_run/gapsBetweenScaffolds -out test_run/gapsBetweenScaffolds/results_MTGLink_DBG_withRef > LOGMTGLinkDBGRunBetweenRef 2> ERRLOGMTGLinkDBGRunBetweenRef  
# if (awk '/\t\* 68-1+:68-2+\n\t\t\* fwd\t2000 bp\tA\n\t\t\* rev\t2000 bp\tA\n\t\* 26939-1+:26939-2+\n\t\t\* fwd\t2004 bp\tB\n\t\t\* rev\t2004 bp\tB/' LOGMTGLinkDBGRunBetweenRef) && (awk '/The gap 27292-1+_27292-2+ was not successfully gap-filled/' LOGMTGLinkDBGRunBetweenRef)
# then
#     echo "MTG-Link DBG w/ ref: Pass"
# else
#     echo "MTG-Link DBG w/ ref: Fail"
# fi

# rm LOGMTGLinkDBGRunBetweenRef
# rm ERRLOGMTGLinkDBGRunBetweenRef
# rm -rf test_run/gapsBetweenScaffolds/results_MTGLink_DBG_withRef

# ### MTG-Link IRO with reference
# /usr/bin/time -v ../mtglink.py IRO -gfa test_run/gapsBetweenScaffolds/test.gfa -bam test_run/gapsBetweenScaffolds/test.bam -fastq test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz -index test_run/gapsBetweenScaffolds/barcoded.bci -s 20 -t 3 -refDir test_run/gapsBetweenScaffolds -out test_run/gapsBetweenScaffolds/results_MTGLink_IRO_withRef > LOGMTGLinkIRORunBetweenRef 2> ERRLOGMTGLinkIRORunBetweenRef  
# if (awk '/\t\* 68-1+:68-2+\n\t\t\* fwd\t1998 bp\tB\n\t\* 26939-1+:26939-2+\n\t\t\* fwd\t2014 bp\tB/' LOGMTGLinkIRORunBetweenRef) && (awk '/The gap 27292-1+_27292-2+ was not successfully gap-filled/' LOGMTGLinkIRORunBetweenRef)
# then
#     echo "MTG-Link IRO w/ ref: Pass"
# else
#     echo "MTG-Link IRO w/ ref: Fail"
# fi

# rm LOGMTGLinkIRORunBetweenRef
# rm ERRLOGMTGLinkIRORunBetweenRef
# rm -rf test_run/gapsBetweenScaffolds/results_MTGLink_IRO_withRef



## Test4: Run of MTG-Link DBG for gaps INTO scaffolds with '--force' parameter

print("TEST4: RUN DBG --force: GAPS INTO SCAFFOLDS")

### MTG-Link DBG --force without reference 
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "--force", "-t", "3", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force"]
commandOut = "LOGMTGLinkDBGRunIntoForce"
commandLog = "ERRLOGMTGLinkDBGRunIntoForce"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k51\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA\n\t\t* fwd_1/2.k41\t1994 bp\tAA\n\t\t* fwd_2/2.k41\t9938 bp\tAA\n\t\t* rev_1/1.k41\t1994 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link DBG --force w/o ref: Pass")
    else:
        print("MTG-Link DBG --force w/o ref: Fail")

subprocess.run("rm LOGMTGLinkDBGRunIntoForce", shell=True)
subprocess.run("rm ERRLOGMTGLinkDBGRunIntoForce", shell=True)
subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force", shell=True)



## Test5: Run of MTG-Link DBG for gaps INTO scaffolds with '--force' parameter and maxLength/minLength

print("TEST5: RUN DBG --force -l -m: GAPS INTO SCAFFOLDS")

### MTG-Link DBG --force -l -m without reference 
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "--force", "-l", "3000", "-m", "2000", "-t", "3", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force_lm"]
commandOut = "LOGMTGLinkDBGRunIntoForce_lm"
commandLog = "ERRLOGMTGLinkDBGRunIntoForce_lm"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k51\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        print("MTG-Link DBG --force -l -m w/o ref: Pass")
    else:
        print("MTG-Link DBG --force -l -m w/o ref: Fail")

subprocess.run("rm LOGMTGLinkDBGRunIntoForce_lm", shell=True)
subprocess.run("rm ERRLOGMTGLinkDBGRunIntoForce_lm", shell=True)
subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force_lm", shell=True)

