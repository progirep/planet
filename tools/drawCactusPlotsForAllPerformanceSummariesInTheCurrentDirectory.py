#!/usr/bin/env python
import os, glob, sys, math

# Settings
# 
PAINTSTYLES = ["color=blue,thick","color=black!50!white,thick","color=red,dashed,thick","color=green!70!black,dotted,thick","color=yellow!60!black","color=red!50!blue","color=green!50!blue","color=green!50!red,thick","color=brown,thick"]
WIDTH = 17
HEIGHT = 8
MAX_TIME = 3600
OVERSHOOT=0.2
NOF_FILES_TICK = 50
MIN_TIME = 0.1
TIME_AXIS_MARKS = [0.1,1,5,10,100,600,1200,3600]

# Time axis computation
def timeToXCoordinate(time):
    time = max(time,MIN_TIME)
    return WIDTH*(math.log(time)-math.log(MIN_TIME)) / (math.log(MAX_TIME)-math.log(MIN_TIME))

# List Performance Files
performanceFiles = glob.glob("*")
performanceFiles = [a for a in performanceFiles if a.find(".")==-1]
performanceFiles.sort()
totalNofFiles = len(glob.glob("../*.rlv"))

# Read Performance Files
performances = {}
for inFile in performanceFiles:

    times = []

    with open(inFile,"r") as inFileReader:
        for line in inFileReader.readlines():
            if line.find("user")!=-1:
                time = float(line[0:line.find("user")])
                times.append(time)
                
    performances[inFile] = times
    
# Write cactus plot   
with open("summary.tex","w") as outFile:    
    print >>outFile,"\\documentclass[a4paper]{scrartcl}"
    print >>outFile,"\\usepackage{tikz}"
    print >>outFile,"\\usepackage[left=1cm,right=1cm,top=2cm,bottom=2cm]{geometry}"
    print >>outFile,"\\title{Version comparison of reluv2}"
    print >>outFile,"\\begin{document}"
    print >>outFile,"\\maketitle"

    # Print legend
    print >>outFile,"\\subsection*{Legend}"
    print >>outFile,"\\begin{tikzpicture}"
    for i, a in enumerate(performances):
        print >>outFile, "\\draw["+PAINTSTYLES[i]+"] (0,",i*0.8,") -- +(2,0) node[anchor=west,color=black] {"+a+"};"
    print >>outFile,"\\end{tikzpicture}"        
    
    # Print cactus plot
    print >>outFile,"\\subsection*{Cactus plot}"
    print >>outFile,"\\centering\\begin{tikzpicture}[scale=0.4]"
    # -> Boundary lines    
    print >>outFile,"\\draw[very thick,->] (0,0) -- (",WIDTH+OVERSHOOT,",0) node[anchor=west] {Time};"
    print >>outFile,"\\draw[very thick,->] (0,0) -- (0,",HEIGHT+OVERSHOOT,") node[anchor=south] {Files};"
    print >>outFile,"\\draw[dashed] (",WIDTH,",0) -- +(0,",HEIGHT,");"
    print >>outFile,"\\draw[dashed] (0,",HEIGHT,") -- +(",WIDTH,",0);"
    # -> Nof files tick
    for i in range(NOF_FILES_TICK,totalNofFiles+1,NOF_FILES_TICK):
        pos = (HEIGHT/float(totalNofFiles))*i
        print >>outFile,"\\draw[semithick] (0.1,",pos,") -- +(-0.2,0) node[anchor=east] {\small "+str(i)+"};"
    # -> Time tick
    for i in TIME_AXIS_MARKS:
        pos = timeToXCoordinate(i)
        print >>outFile,"\\draw[semithick] (",pos,",0.1) -- +(0,-0.2) node[anchor=north] {\small "+str(i)+"};"
    # -> Cactus plot
    for i,filename in enumerate(performances):
        allTimes = performances[filename]
        style = PAINTSTYLES[i]
        allTimes.sort()
        print >>outFile,"\\draw["+style+"] (0,0)",
        nofFilesSoFar = 0
        for t in allTimes:
            nofFilesSoFar += 1
            xPos = timeToXCoordinate(t)
            yPos = (HEIGHT/float(totalNofFiles))*nofFilesSoFar
            print >>outFile,"-| (",xPos,",",yPos,")",
        print >>outFile,";"
    print >>outFile,"\\end{tikzpicture}\n\n"        
    
    
    
    
    
    
    
    print >>outFile,"\\end{document}"
    
                
    


