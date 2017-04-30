#!/usr/bin/python2
import os, sys

for inFile in sys.argv[1:]:

    sumTimes = 0.0
    with open(inFile,"r") as inFileReader:
        for line in inFileReader.readlines():
            if line.find("user")!=-1:
                time = float(line[0:line.find("user")])
                sumTimes += time
    
    print str(sumTimes)+"user"
