#!/usr/bin/env python2
import os, glob, sys, subprocess
import PIL.Image
from joblib import Parallel, delayed

if len(sys.argv)<3:
    print >>sys.stderr, "Error: Wanted an RLV filename (output by json_network_to_rlv_translator) and a path with the testing PNG images"
    sys.exit(1)
    
rlvFile = sys.argv[1]
testingDir = sys.argv[2]

# Read RLV Lines
with open(rlvFile,"r") as inRLV:
    rlvLines = inRLV.readlines()
    
# Iterate over digits
correct = 0
total = 0

for digit in range(0,10):
    print "\nDigit: ",digit
    allPngFiles = glob.glob(testingDir+"/"+str(digit)+"/*.png")
    
    def workOnPNG(pngFile):
        im = PIL.Image.open(pngFile) #Can be many different formats.
        pix = im.load()
        assert im.size==(28,28)
        pixels = [pix[x,y] for y in range(0,28) for x in range(0,28)]
        
        # Open RELU-Verifier
        verifierProcess = subprocess.Popen("../ReLUVerifier/src2/reluv2 /dev/stdin",shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        # outtaFile = open("/tmp/iniCheck.rlv","w")
        for a in rlvLines:
            verifierProcess.stdin.write(a)
            # outtaFile.write(a)
        for i in range(0,len(pixels)):
            pixelValue = pixels[i]/256.0
            verifierProcess.stdin.write("Assert <= "+str(pixelValue-0.0001)+" "+"1.0 inX"+str(i)+"\n")
            verifierProcess.stdin.write("Assert >= "+str(pixelValue+0.0001)+" "+"1.0 inX"+str(i)+"\n")
            # outtaFile.write("Assert <= "+str(pixelValue)+" "+"1.0 inX"+str(i)+"\n")
            # outtaFile.write("Assert >= "+str(pixelValue)+" "+"1.0 inX"+str(i)+"\n")
        verifierProcess.stdin.close()
        # outtaFile.close()
        
        foundSATLine = False
        foundValuationLine = False
        values = {}
        for a in verifierProcess.stdout.readlines():
            a = a.strip()
            if a=="SAT":
                foundSATLine = True
            elif a=="Valuation:":
                foundValuationLine = True
            elif a.startswith("- ") and foundValuationLine:
                parts = a.split(" ")
                assert parts[0]=="-"
                assert parts[3]=="/"
                values[parts[1][0:len(parts[1])-1]] = float(parts[2])

        assert verifierProcess.wait()==0
        assert foundSATLine==True
        
        bestDigit = -1
        bestValue = None
        for i in range(0,10):
            if (bestValue==None) or values["outX"+str(i)]>bestValue:
                bestValue = values["outX"+str(i)]
                bestDigit = i
        
        sys.stdout.write(str(bestDigit)+",")
        sys.stdout.flush()
        return bestDigit
    
    allResults = Parallel(n_jobs=4, backend="threading")(delayed(workOnPNG)(a) for a in allPngFiles)
    total += len(allResults)
    for a in allResults:
        if a==digit:
            correct += 1
            
print "Correct were ", correct, "out of", total
