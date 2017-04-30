#!/usr/bin/env python2
import os, glob, sys, subprocess
import PIL.Image
from joblib import Parallel, delayed

if len(sys.argv)<3:
    print >>sys.stderr, "Error: Wanted an RLV filename, a mode, and parameters for the modes"
    print >>sys.stderr, "Modes are:"
    print >>sys.stderr, "- GIVE <digit> - Obtain an arbitrary image that yields the number"
    print >>sys.stderr, "- GIVESTRONG <digit> <minDifference>- Obtain an arbitrary image that strongly yields the number, i.e., the other digit values in the output layer must be at least <minDifference> smaller."
    print >>sys.stderr, "- ROBUST <digitFile> <targetDigit> <maxDifferencePerPixel> <maxUnsmoothnessInNoise>- Obtain a digit image that is close to a given one, that resolves to the given target digit, where every pixel is at most maxDifferencePerPixel away from the initial image, and the maximal noise difference between two adjacent pixels is maxUnsmoothnessInNoise. The last two parameters should be >=0 and <=1 (such as, e.g., 0.05 for 5% deviation)"
    sys.exit(1)
    
rlvFile = sys.argv[1]
mode = sys.argv[2]

# Read RLV Lines
with open(rlvFile,"r") as inRLV:
    rlvLines = inRLV.readlines()


if mode=="ANY":

    # Set the boundaries
    verifierProcess = subprocess.Popen("../../src/planet /dev/stdin",shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    for a in rlvLines:
        verifierProcess.stdin.write(a)
    for i in range(0,28*28):
        verifierProcess.stdin.write("Assert <= 0.0 1.0 inX"+str(i)+"\n")
        verifierProcess.stdin.write("Assert >= 0.0 1.0 inX"+str(i)+"\n")
        
    verifierProcess.stdin.close()

    foundSATLine = False
    foundValuationLine = False
    values = {}
    for a in verifierProcess.stdout.readlines():
        sys.stdout.write(a)
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

    if not foundSATLine:
        print "No digit found."
    else:
        outFileName = "/tmp/outImage"+str(os.getpid())+"-any.png"
        outImage = PIL.Image.new("L", (28, 28))
        for y in range(0,28):
            for x in range(0,28):    
                outImage.putpixel((x,y),int(256*values["inX"+str(y*28+x)]))
                
        outImage.save(outFileName)
        print "Result is found in image file:",outFileName


# Give mode
elif mode=="GIVE" or mode=="GIVESTRONG":

    if len(sys.argv)<4:
        print >>sys.stderr, "Error: GIVE and GIVESTRING mode requires a digit number."
        sys.exit(1)
    digit = int(sys.argv[3])
    
    if mode=="GIVESTRONG":
        if len(sys.argv)<5:
            print >>sys.stderr, "Error: GIVESTRING requires a minimum difference."
            sys.exit(1)
        minDifference = -1*float(sys.argv[4])
    else:        
        minDifference = -0.0000001    
    

    # Set the boundaries
    verifierProcess = subprocess.Popen("../../src/planet /dev/stdin",shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    for a in rlvLines:
        verifierProcess.stdin.write(a)
    for i in range(0,28*28):
        verifierProcess.stdin.write("Assert <= 0.0 1.0 inX"+str(i)+"\n")
        verifierProcess.stdin.write("Assert >= 1.0 1.0 inX"+str(i)+"\n")
        
    # Set the output
    for i in range(0,10):
        if i!=digit:
            verifierProcess.stdin.write("Assert >= "+str(minDifference)+" 1.0 outX"+str(i)+" -1.0 outX"+str(digit)+"\n")
    verifierProcess.stdin.close()

    foundSATLine = False
    foundValuationLine = False
    values = {}
    for a in verifierProcess.stdout.readlines():
        sys.stdout.write(a)
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

    if not foundSATLine:
        print "No digit found."
    else:
        outFileName = "/tmp/outImage"+str(os.getpid())+"-"+str(digit)+str(minDifference)+".png"
        outImage = PIL.Image.new("L", (28, 28))
        for y in range(0,28):
            for x in range(0,28):    
                outImage.putpixel((x,y),int(256*values["inX"+str(y*28+x)]))
                
        outImage.save(outFileName)
        print "Result is found in image file:",outFileName


# ROBUST mode
elif mode=="ROBUST":

    if len(sys.argv)<7:
        print >>sys.stderr, "Error: ROBUST needs many parameters."
        sys.exit(1)
    digitFile = sys.argv[3]
    targetDigit = int(sys.argv[4])
    maxDifferencePerPixel = float(sys.argv[5])
    maxUnsmoothnessInNoise = float(sys.argv[6])
    
    # Read image file
    im = PIL.Image.open(digitFile) #Can be many different formats.
    pix = im.load()
    assert im.size==(28,28)
    pixels = [pix[x,y] for y in range(0,28) for x in range(0,28)]
    
    # Set the boundaries
    verifierProcess = subprocess.Popen("../../src/planet /dev/stdin",shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE)
    for a in rlvLines:
        verifierProcess.stdin.write(a)
    for i in range(0,28*28):
        # Outermost pixels need to be excluded
        x = i % 28
        y = int(i / 28)
        if x<3 or x>24 or y<3 or y>24:
            border = 0.0
        else:
            border = maxDifferencePerPixel
        mini = max(0.0,pixels[i]/256.0-border)
        maxi = min(1.0,pixels[i]/256.0+border)
        verifierProcess.stdin.write("Assert <= "+str(mini)+" 1.0 inX"+str(i)+"\n")
        verifierProcess.stdin.write("Assert >= "+str(maxi)+" 1.0 inX"+str(i)+"\n")
        
    # Set the output
    for i in range(0,10):
        if i!=targetDigit:
            verifierProcess.stdin.write("Assert >= -0.000001 1.0 outX"+str(i)+" -1.0 outX"+str(targetDigit)+"\n")
            
    # Set the smoothness
    if maxUnsmoothnessInNoise<1.0:
        for x in range(0,28):
            for y in range(0,28):
                # Smooth down
                if (y<27):
                    pixelDiff = (pixels[y*28+x]-pixels[(y+1)*28+x])/256.0
                    verifierProcess.stdin.write("Assert <= "+str(pixelDiff-maxUnsmoothnessInNoise)+" 1.0 inX"+str(y*28+x)+" -1.0 inX"+str((y+1)*28+x)+"\n")
                    verifierProcess.stdin.write("Assert >= "+str(pixelDiff+maxUnsmoothnessInNoise)+" 1.0 inX"+str(y*28+x)+" -1.0 inX"+str((y+1)*28+x)+"\n")
                # Smooth right
                if (x<27):
                    pixelDiff = (pixels[y*28+x]-pixels[y*28+x+1])/256.0
                    verifierProcess.stdin.write("Assert <= "+str(pixelDiff-maxUnsmoothnessInNoise)+" 1.0 inX"+str(y*28+x)+" -1.0 inX"+str(y*28+x+1)+"\n")
                    verifierProcess.stdin.write("Assert >= "+str(pixelDiff+maxUnsmoothnessInNoise)+" 1.0 inX"+str(y*28+x)+" -1.0 inX"+str(y*28+x+1)+"\n")
                
    # Done with the input instance    
    verifierProcess.stdin.close()

    foundSATLine = False
    foundValuationLine = False
    values = {}
    for a in verifierProcess.stdout.readlines():
        sys.stdout.write(a)
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

    if not foundSATLine:
        print "No digit found."
    else:
        outFileName = "/tmp/outImage"+str(os.getpid())+"-"+str(targetDigit)+"-"+str(maxDifferencePerPixel)+"-"+str(maxUnsmoothnessInNoise)+".png"
        outImage = PIL.Image.new("L", (28, 28))
        for y in range(0,28):
            for x in range(0,28):    
                outImage.putpixel((x,y),int(256*values["inX"+str(y*28+x)]))
                
        outImage.save(outFileName)
        print "Result is found in image file:",outFileName



else:
    # Unknown mode.
    print >>sys.stderr, "Unknown 'prodNetwork' operation mode: ",mode
    sys.exit(1)
