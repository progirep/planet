#!/usr/bin/python2
#
# Build a makefile for checking many RLV cases in parallel
import os, sys, glob
    
if len(sys.argv)!=2:
    print >>sys.stderr, "Error: Too many or too few parameters passed. Expecting only the parameter to PLANET (--smtlibApprox or --smtlib)"
    sys.exit(1)
    
scriptpath = os.path.dirname(os.path.realpath(__file__))

allRLVFiles = list(glob.glob("*.rlv"))
allRLVFiles.sort()
option = sys.argv[1]

with open("Makefile","w") as outFile:
    outFile.write("default: "+" ".join([".out/"+a+".ytxt" for a in allRLVFiles]))
    outFile.write("\n\t@mkdir -p performanceComparisons\n")
    outFile.write("\t@cat .out/*.ytxt > \"performanceComparisons/YICES`(git show | grep Date)`\"\n\n")

    for a in allRLVFiles:
        outFile.write(".out/"+a+".ytxt:\n")
        outFile.write("\t@mkdir -p .out\n")
        outFile.write("\t@echo \"Processing: "+a+"\"\n")
        outFile.write("\ttime "+scriptpath+"/../src/planet "+option+" "+a+" > .out/"+a+".smt2 2> .out/"+a+".ytxt \n")
        outFile.write("\t@echo \"Processing Part 2: "+a+"\"\n")
        if sys.argv[1].find("pprox")>0:
            outFile.write("\t@(time timeout 3600 $$YICES/bin/yices-smt2 .out/"+a+".smt2 >> .out/"+a+".ytxt 2>&1) || if [ $$? -eq 124 ]; then echo \"--> "+a+" timed out\"; else false; fi\n")
        else:
            outFile.write("\t@(time timeout 3600 $$YICES/bin/yices-smt2 .out/"+a+".smt2 > .out/"+a+".ytxt 2>&1) || if [ $$? -eq 124 ]; then echo \"--> "+a+" timed out\"; else false; fi\n")
        outFile.write("\tcat .out/"+a+".ytxt | grep \"^Command exited with non-zero status 124\\|^")
        if a.endswith("_SAT.rlv"):
            outFile.write("sat\"\n")
        elif a.endswith("_UNSAT.rlv"):
            outFile.write("unsat\"\n")
        else:
            raise "File does not end with _SAT or _UNSAT: "+a
        outFile.write("\n\n")
        
