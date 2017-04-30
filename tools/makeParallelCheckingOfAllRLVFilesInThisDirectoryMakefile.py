#!/usr/bin/python2
#
# Build a makefile for checking many RLV cases in parallel
import os, sys, glob

if len(sys.argv)<2:
    print >>sys.stderr, "Error: Expected the path to the verifier as parameter."
    sys.exit(1)
    
if len(sys.argv)>3:
    print >>sys.stderr, "Error: Too many parameters passed. Expecting the revifier command and optionally a file to check."
    sys.exit(1)
    
verifier = sys.argv[1]
scriptpath = os.path.dirname(os.path.realpath(__file__))

allRLVFiles = list(glob.glob("*.rlv"))
allRLVFiles.sort()

with open("Makefile","w") as outFile:
    outFile.write("default: "+" ".join([".out/"+a+".txt" for a in allRLVFiles]))
    outFile.write("\n\t@mkdir -p performanceComparisons\n")
    outFile.write("\t@cat .out/*.txt > \"performanceComparisons/`(git show | grep Date)`\"\n\n")

    for a in allRLVFiles:
        outFile.write(".out/"+a+".txt:\n")
        outFile.write("\t@mkdir -p .out\n")
        outFile.write("\t@echo \"Processing: "+a+"\"\n")
        outFile.write("\t@(time timeout 3600 "+scriptpath+"/"+"testBenchmarksInThisDirectory.py "+verifier+" "+a+" > .out/"+a+".txt 2>&1) || if [ $$? -eq 124 ]; then echo \"--> "+a+" timed out\"; else false; fi \n\n")
        
