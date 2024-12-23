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
        passed = True
        print("MTG-Link DBG: Pass")
    else:
        passed = False
        print("MTG-Link DBG: Fail")

if passed:
    subprocess.run("rm LOGMTGLinkDBGInstall", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGInstall", shell=True)
    subprocess.run("rm -rf test_install/results_MTGLink_DBG", shell=True)


# ### MTG-Link IRO
# command = ["/usr/bin/time", "-v", "../mtglink.py", "IRO", \
#             "-gfa", "test_install/test.gfa", "-bam", "test_install/test.bam", \
#             "-fastq", "test_install/reads.sorted.fastq.gz", "-index", "test_install/barcoded.bci", \
#             "-out", "test_install/results_MTGLink_IRO"]
# commandOut = "LOGMTGLinkIROInstall"
# commandLog = "ERRLOGMTGLinkIROInstall"
# with open(commandOut, "a") as out, open(commandLog, "a") as log:
#     subprocess.run(command, stdout=out, stderr=log)

# expectedString = "In total, 1 targets were successfully assembled:\n\t* 8_0-3152098-L+:8_3153098-6305195-R+\n\t\t* fwd.k0\t2000 bp\tAA"
# with open(commandOut) as out:
#     if expectedString in out.read():
#         print("MTG-Link IRO: Pass")
#     else:
#         print("MTG-Link IRO: Fail")

# subprocess.run("rm LOGMTGLinkIROInstall", shell=True)
# subprocess.run("rm ERRLOGMTGLinkIROInstall", shell=True)
# subprocess.run("rm -rf test_install/results_MTGLink_IRO", shell=True)



## Test2: Run of MTG-Link for gaps INTO scaffolds

print("TEST2: RUN: GAPS INTO SCAFFOLDS")

### MTG-Link DBG
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "31", "21", "-t", "4", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG"]
commandOut = "LOGMTGLinkDBGRunInto"
commandLog = "ERRLOGMTGLinkDBGRunInto"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        passed = True
        print("MTG-Link DBG: Pass")
    else:
        passed = False
        print("MTG-Link DBG: Fail")

# # We need these files for the test with the '-bxuDir' arg, so do not remove it now. 
# subprocess.run("rm LOGMTGLinkDBGRunInto", shell=True)
# subprocess.run("rm ERRLOGMTGLinkDBGRunInto", shell=True)
# subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG", shell=True)


# ### MTG-Link IRO
# command = ["/usr/bin/time", "-v", "../mtglink.py", "IRO", \
#             "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
#             "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
#             "-s", "20", "-t", "3", \
#             "-out", "test_run/gapsIntoScaffolds/results_MTGLink_IRO"]
# commandOut = "LOGMTGLinkIRORunInto"
# commandLog = "ERRLOGMTGLinkIRORunInto"
# with open(commandOut, "a") as out, open(commandLog, "a") as log:
#     subprocess.run(command, stdout=out, stderr=log)

# expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd.k0\t1998 bp\tAB\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd.k0\t2014 bp\tAA"
# with open(commandOut) as out:
#     if expectedString in out.read():
#         print("MTG-Link IRO: Pass")
#     else:
#         print("MTG-Link IRO: Fail")

# subprocess.run("rm LOGMTGLinkIRORunInto", shell=True)
# subprocess.run("rm ERRLOGMTGLinkIRORunInto", shell=True)
# subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_IRO", shell=True)




## Test3: Run of MTG-Link for gaps BETWEEN scaffolds

print("TEST3: RUN: GAPS BETWEEN SCAFFOLDS")

### MTG-Link DBG
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsBetweenScaffolds/test.gfa", "-bam", "test_run/gapsBetweenScaffolds/test.bam", \
            "-fastq", "test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsBetweenScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "31", "21", "-t", "4", \
            "-out", "test_run/gapsBetweenScaffolds/results_MTGLink_DBG"]
commandOut = "LOGMTGLinkDBGRunBetween"
commandLog = "ERRLOGMTGLinkDBGRunBetween"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68-1+:68-2+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t* 26939-1+:26939-2+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        passed = True
        print("MTG-Link DBG: Pass")
    else:
        passed = False
        print("MTG-Link DBG: Fail")

if passed:
    subprocess.run("rm LOGMTGLinkDBGRunBetween", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGRunBetween", shell=True)
    subprocess.run("rm -rf test_run/gapsBetweenScaffolds/results_MTGLink_DBG", shell=True)


# ### MTG-Link IRO
# command = ["/usr/bin/time", "-v", "../mtglink.py", "IRO", \
#             "-gfa", "test_run/gapsBetweenScaffolds/test.gfa", "-bam", "test_run/gapsBetweenScaffolds/test.bam", \
#             "-fastq", "test_run/gapsBetweenScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsBetweenScaffolds/barcoded.bci", \
#             "-s", "20", "-t", "3", \
#             "-out", "test_run/gapsBetweenScaffolds/results_MTGLink_IRO"]
# commandOut = "LOGMTGLinkIRORunBetween"
# commandLog = "ERRLOGMTGLinkIRORunBetween"
# with open(commandOut, "a") as out, open(commandLog, "a") as log:
#     subprocess.run(command, stdout=out, stderr=log)

# expectedString = "In total, 2 targets were successfully assembled:\n\t* 68-1+:68-2+\n\t\t* fwd.k0\t1998 bp\tAB\n\t* 26939-1+:26939-2+\n\t\t* fwd.k0\t2014 bp\tAA"
# with open(commandOut) as out:
#     if expectedString in out.read():
#         print("MTG-Link IRO: Pass")
#     else:
#         print("MTG-Link IRO: Fail")

# subprocess.run("rm LOGMTGLinkIRORunBetween", shell=True)
# subprocess.run("rm ERRLOGMTGLinkIRORunBetween", shell=True)
# subprocess.run("rm -rf test_run/gapsBetweenScaffolds/results_MTGLink_IRO", shell=True)




## Test4: Run of MTG-Link DBG for gaps INTO scaffolds with '-bxuDir' parameter

print("TEST4: RUN DBG -bxuDir: GAPS INTO SCAFFOLDS")

### MTG-Link DBG -bxuDir
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-bxuDir", "test_run/gapsIntoScaffolds/results_MTGLink_DBG/read_subsampling/", \
            "-k", "61", "51", "41", "31", "21", "-t", "4", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG_bxuDir"]
commandOut = "LOGMTGLinkDBGRunIntoBXUDir"
commandLog = "ERRLOGMTGLinkDBGRunIntoBXUDir"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 2 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        passed = True
        print("MTG-Link DBG -bxuDir: Pass")
    else:
        passed = False
        print("MTG-Link DBG -bxuDir: Fail")

if passed:
    subprocess.run("rm LOGMTGLinkDBGRunInto", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGRunInto", shell=True)
    subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG", shell=True)
    subprocess.run("rm LOGMTGLinkDBGRunIntoBXUDir", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGRunIntoBXUDir", shell=True)
    subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_bxuDir", shell=True)



## Test5: Run of MTG-Link DBG for gaps INTO scaffolds with '--multiple' parameter

print("TEST5: RUN DBG --multiple: GAPS INTO SCAFFOLDS")

### MTG-Link DBG --multiple
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "31", "21", "--multiple", "-t", "4", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG_Multiple"]
commandOut = "LOGMTGLinkDBGRunIntoMultiple"
commandLog = "ERRLOGMTGLinkDBGRunIntoMultiple"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 3 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA\n\t* 27358_0-171884-L+:27358_172884-344768-R+\n\t\t* fwd_1/2.k61\t1756 bp\tAA\n\t\t* fwd_2/2.k61\t2000 bp\tAA\n\t\t* rev_1/2.k61\t2000 bp\tAA\n\t\t* rev_2/2.k61\t1756 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        passed = True
        print("MTG-Link DBG --multiple: Pass")
    else:
        passed = False
        print("MTG-Link DBG --multiple: Fail")

if passed:
    subprocess.run("rm LOGMTGLinkDBGRunIntoMultiple", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGRunIntoMultiple", shell=True)
    subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_Multiple", shell=True)



## Test6: Run of MTG-Link DBG for gaps INTO scaffolds with '--force' parameter

print("TEST6: RUN DBG --force: GAPS INTO SCAFFOLDS")

### MTG-Link DBG --force
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "--force", "-t", "4", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force"]
commandOut = "LOGMTGLinkDBGRunIntoForce"
commandLog = "ERRLOGMTGLinkDBGRunIntoForce"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 3 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k51\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA\n\t* 27358_0-171884-L+:27358_172884-344768-R+\n\t\t* fwd_1/2.k61\t1756 bp\tAA\n\t\t* fwd_2/2.k61\t2000 bp\tAA\n\t\t* rev_1/2.k61\t2000 bp\tAA\n\t\t* rev_2/2.k61\t1756 bp\tAA\n\t\t* fwd_1/2.k51\t1756 bp\tAA\n\t\t* fwd_2/2.k51\t2000 bp\tAA\n\t\t* rev_1/2.k51\t2000 bp\tAA\n\t\t* rev_2/2.k51\t1756 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        passed = True
        print("MTG-Link DBG --force: Pass")
    else:
        passed = False
        print("MTG-Link DBG --force: Fail")

if passed:
    subprocess.run("rm LOGMTGLinkDBGRunIntoForce", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGRunIntoForce", shell=True)
    subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force", shell=True)



## Test7: Run of MTG-Link DBG for gaps INTO scaffolds with '--force' parameter and maxLength/minLength

print("TEST7: RUN DBG --force -l -m: GAPS INTO SCAFFOLDS")

### MTG-Link DBG --force -l -m
command = ["/usr/bin/time", "-v", "../mtglink.py", "DBG", \
            "-gfa", "test_run/gapsIntoScaffolds/test.gfa", "-bam", "test_run/gapsIntoScaffolds/test.bam", \
            "-fastq", "test_run/gapsIntoScaffolds/reads.sorted.fastq.gz", "-index", "test_run/gapsIntoScaffolds/barcoded.bci", \
            "-k", "61", "51", "41", "--force", "-l", "3000", "-m", "2000", "-t", "4", \
            "-out", "test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force_lm"]
commandOut = "LOGMTGLinkDBGRunIntoForce_lm"
commandLog = "ERRLOGMTGLinkDBGRunIntoForce_lm"
with open(commandOut, "a") as out, open(commandLog, "a") as log:
    subprocess.run(command, stdout=out, stderr=log)

expectedString = "In total, 3 targets were successfully assembled:\n\t* 68_0-222132-L+:68_223132-445265-R+\n\t\t* fwd_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k61\t2000 bp\tAA\n\t\t* rev_1/1.k51\t2000 bp\tAA\n\t* 26939_0-464805-L+:26939_465805-930610-R+\n\t\t* fwd_1/1.k51\t2004 bp\tAA\n\t\t* rev_1/1.k51\t2004 bp\tAA\n\t* 27358_0-171884-L+:27358_172884-344768-R+\n\t\t* fwd_2/2.k61\t2000 bp\tAA\n\t\t* rev_1/2.k61\t2000 bp\tAA\n\t\t* fwd_2/2.k51\t2000 bp\tAA\n\t\t* rev_1/2.k51\t2000 bp\tAA\n\t\t* rev_2/2.k41\t2000 bp\tAA"
with open(commandOut) as out:
    if expectedString in out.read():
        passed = True
        print("MTG-Link DBG --force -l -m: Pass")
    else:
        passed = False
        print("MTG-Link DBG --force -l -m: Fail")

if passed:
    subprocess.run("rm LOGMTGLinkDBGRunIntoForce_lm", shell=True)
    subprocess.run("rm ERRLOGMTGLinkDBGRunIntoForce_lm", shell=True)
    subprocess.run("rm -rf test_run/gapsIntoScaffolds/results_MTGLink_DBG_Force_lm", shell=True)

