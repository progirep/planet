#!/usr/bin/python2
import os, sys

for inFile in sys.argv[1:]:

    sumTimes = 0.0
    foundAbortString = False
    with open(inFile,"r") as inFileReader:
        for line in inFileReader.readlines():
            if line.find("user")!=-1:
                time = float(line[0:line.find("user")])
                sumTimes += time
            if line.find("Command exited with non-zero status 124")!=-1:
                foundAbortString = True
    
    if not foundAbortString:
        print str(sumTimes)+"user"
