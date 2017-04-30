#!/usr/bin/env python2
import os, glob, sys, subprocess
import PIL.Image
from joblib import Parallel, delayed

if len(sys.argv)<2:
    print >>sys.stderr, "Error: Wanted an PNG filename"
    sys.exit(1)
    
pngFile = sys.argv[1]
inImage = PIL.Image.open(pngFile)

print "\\begin{tikzpicture}[xscale=0.09,yscale=-0.09]"
xSize = inImage.width
ySize = inImage.height
for y in range(0,ySize):
    for x in range(0,xSize):
        print "\\path[fill=black!"+str(inImage.getpixel((x,y))/2.550)+"!white] (",x,",",y,") rectangle +(1,1);"
#for y in range(0,ySize):
#    print "\\draw[line width=0.03pt,color=black!10!white] (0,",y,") -- +(",xSize,",0);"
#for x in range(0,xSize):
#    print "\\draw[line width=0.03pt,color=black!10!white] (",x,",0) -- +(0,",ySize,");"

print "\\draw[semithick] (0,0) rectangle (",xSize,",",ySize,");"
print "\\end{tikzpicture}"

