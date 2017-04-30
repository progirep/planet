#!/usr/bin/env python2
#
# This program takes a JSON file that describes a network for classification that has been learned from a CSV file
# and the CSV file itself. It then checks if the Network has an accurracy of 1 (i.e., whether every data point is classified correctly)
# and what the robustness against deviations in every input dimension is.
#
# It uses other tools from the "tools" directory to achieve this, including the RLV compressor.
import os, sys,subprocess, time

VERIFIER = "/../../src/planet"
ABORT_ERROR_MARGIN = 0.002
START_ERROR_MARGIN = 0.05

if len(sys.argv)<3:
    print >>sys.stderr, "Error: Need two file names as parameters: the JSON-File and the CSV file."
    sys.exit(1)

# Read CSV file    
with open(sys.argv[2],"r") as inFile:
    csvHeader = inFile.readline()
    csvLines = [a.strip() for a in inFile.readlines()]
    
# Process CSV file
nofClasses = 0
for i in range(0,len(csvLines)):
    parts = csvLines[i].split(",")
    parts[-1] = int(parts[-1])
    nofClasses = max(nofClasses,parts[-1])
    for j in range(0,len(parts)-1):
        parts[j] = float(parts[j])
    csvLines[i] = parts
    
    
# Translate to RLV without main constraint
tmpFilenameUnpreprocessedRLVFile = "/tmp/RLVUnpreprocessed-"+str(os.getpid())+".rlv"
tmpFilenamePreprocessedRLVFile = "/tmp/RLVProcessed-"+str(os.getpid())+".rlv"
scriptPath = os.path.dirname(sys.argv[0])
translatorToRLVProcess = subprocess.Popen([scriptPath+"/../../tools/json_network_to_rlv_translator.py",sys.argv[1]],stdout=subprocess.PIPE)
allRLVLines = [a.strip() for a in translatorToRLVProcess.stdout.readlines()]
assert translatorToRLVProcess.wait()==0

nofCSVLinesProcessed = 0
for csvLine in csvLines:
    nofCSVLinesProcessed += 1
    if nofCSVLinesProcessed>100:
        print >>sys.stderr, "Processed 100 cases."
        sys.exit(0)
    print " ".join([str(a) for a in csvLine]),"->",
    sys.stdout.flush()

    # Build a "precise instance"
    isAnySAT = False
    for j in range(0,nofClasses+1):
        if j!=csvLine[-1]:   
            with open(tmpFilenameUnpreprocessedRLVFile,"w") as outFile:
                for line in allRLVLines:
                    outFile.write(line+"\n")
                for i in range(0,len(csvLine)-1):
                    outFile.write("Assert <= "+str(csvLine[i])+" 1.0 inX"+str(i)+"\n")
                    outFile.write("Assert >= "+str(csvLine[i])+" 1.0 inX"+str(i)+"\n")
                outFile.write("Assert >= 0.0 -1.0 outX"+str(j)+" 1.0 outX"+str(csvLine[-1])+"\n")
            
            # Now Preprocess    
            preprocessorProcess = subprocess.Popen([scriptPath+"/../../tools/rlvCompressor.py",tmpFilenameUnpreprocessedRLVFile],stdout=subprocess.PIPE)
            allPreprocessedLines = preprocessorProcess.stdout.readlines()
            assert preprocessorProcess.wait()==0
            with open(tmpFilenamePreprocessedRLVFile,"w") as outFile:
                for a in allPreprocessedLines:
                    outFile.write(a)
                
            # Now check with the Checker.... 
            modelCheckerProcess = subprocess.Popen([scriptPath+VERIFIER,tmpFilenamePreprocessedRLVFile],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            allMCLines = modelCheckerProcess.stdout.readlines()
            if modelCheckerProcess.wait()!=0:
                print >>sys.stderr, "Error executing verifier:"
                for a in allMCLines:
                    print a
                sys.exit(1)
            isSAT = None
            for line in allMCLines:
                if line.startswith("SAT"):
                    isSAT = True
                elif line.startswith("UNSAT"):
                    isSAT = False
            sys.stdout.flush()
            assert not isSAT is None
            isAnySAT = isAnySAT or isSAT
            
     # Finally print result
    if isAnySAT:
        print "Found Error"
    else:
        print "OK,",
        sys.stdout.flush()
        
        # Now compute robustness
        errorMarginMin = 0.0
        errorMarginMax = START_ERROR_MARGIN
        while ((errorMarginMax-errorMarginMin) > ABORT_ERROR_MARGIN):
        
            errorMargin = (errorMarginMax+errorMarginMin)/2.0
            foundError = False
            for j in range(0,nofClasses+1):
                if j!=csvLine[-1]:   
                    with open(tmpFilenameUnpreprocessedRLVFile,"w") as outFile:
                        for line in allRLVLines:
                            outFile.write(line+"\n")
                        for i in range(0,len(csvLine)-1):
                            outFile.write("Assert <= "+str(csvLine[i]-errorMargin)+" 1.0 inX"+str(i)+"\n")
                            outFile.write("Assert >= "+str(csvLine[i]+errorMargin)+" 1.0 inX"+str(i)+"\n")
                        outFile.write("Assert >= 0.0 -1.0 outX"+str(j)+" 1.0 outX"+str(csvLine[-1])+"\n")


                    # Now Preprocess    
                    preprocessorProcess = subprocess.Popen([scriptPath+"/../../tools/rlvCompressor.py",tmpFilenameUnpreprocessedRLVFile],stdout=subprocess.PIPE)
                    allPreprocessedLines = preprocessorProcess.stdout.readlines()
                    assert preprocessorProcess.wait()==0
                    with open(tmpFilenamePreprocessedRLVFile,"w") as outFile:
                        for a in allPreprocessedLines:
                            outFile.write(a)

                        
                    # Now check with the Checker.... 
                    startingTime = time.time()
                    modelCheckerProcess = subprocess.Popen([scriptPath+VERIFIER,tmpFilenamePreprocessedRLVFile],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                    allMCLines = modelCheckerProcess.stdout.readlines()
                    timeDuration = time.time() - startingTime
                    
                    assert modelCheckerProcess.wait()==0
                    isSAT = None
                    for line in allMCLines:
                        if line.startswith("SAT"):
                            isSAT = True
                        elif line.startswith("UNSAT"):
                            isSAT = False
                    sys.stdout.flush()
                    assert not isSAT is None
                    foundError = foundError or isSAT

                    # Write benchmark file if applicable.
                    # if (timeDuration<6000.0) and (timeDuration>0.05):
                    if True:
                        satPostfix = "SAT"
                        if isSAT==False:
                            satPostfix = "UNSAT"
                        os.system("cp "+tmpFilenamePreprocessedRLVFile+" /tmp/reluBenchmark"+str(timeDuration)+"s_"+satPostfix+".rlv")

        
            if foundError:
                errorMarginMax = errorMargin
            else:
                errorMarginMin = errorMargin
        print "Error Margin >=",errorMarginMin
