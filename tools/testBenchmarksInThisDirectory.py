#!/usr/bin/env python2
import os, glob, sys, subprocess

if len(sys.argv)<2:
    print >>sys.stderr, "Error: Expected the path to the verifier as parameter."
    sys.exit(1)
    
if len(sys.argv)>3:
    print >>sys.stderr, "Error: Too many parameters passed. Expecting the revifier command and optionally a file to check."
    sys.exit(1)
    
verifier = sys.argv[1]

if len(sys.argv)==3:
    allRLVFiles = list(glob.glob(sys.argv[2]))
else:
    allRLVFiles = list(glob.glob("*.rlv"))
allRLVFiles.sort()

# Find RELUV binary
location = sys.argv[0]
while location[-1]!="/" and len(location)>0:
    location = location[0:len(location)-1]

# Iterate over the benchmarks
for benchmark in allRLVFiles:
    assert benchmark.endswith(".rlv")
    sys.stdout.write(benchmark)
    sys.stdout.flush()
    benchmarkMainPart = benchmark[0:len(benchmark)-4]
    
    # Run RLV
    process = subprocess.Popen([verifier,benchmark],stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    resultFound = False
    result = None
    totalErrorLine = None
    otherLines = []
    for line in process.stdout.readlines():
        otherLines.append(line)    
        line = line + "\n"
        if line.startswith("SAT\n"):
            resultFound = True
            result = True
        elif line.startswith("UNSAT\n"):
            resultFound = True
            result = False
        elif line.startswith("Total error:"):
            totalErrorLine = line.strip()
    errorCode = process.wait()
    if errorCode!=0:
        print >>sys.stderr, "Error running "+verifier+" "+benchmark+ " -- Error Code:",errorCode
        print >>sys.stderr, "Output:\n"+"".join(otherLines)
        print >>sys.stderr, "That was an error running "+verifier+" "+benchmark+ " -- Error Code:",errorCode        
        sys.exit(1)
    
    if not resultFound:
        print >>sys.stderr, "Error: No result for "+benchmark+" found."
        sys.exit(1)
    
    if benchmarkMainPart.endswith("_SAT"):
        if not result:
            print >>sys.stderr, "Error: Differing results for "+benchmark+" found."
            sys.exit(1)
    elif benchmarkMainPart.endswith("_UNSAT"):
        if result:
            print >>sys.stderr, "Error: Differing results for "+benchmark+" found."
            sys.exit(1)
    else:
        print "Note: No expected result for "+benchmark+" found."
    
    print " -> OK."
    if totalErrorLine!=None:
        print totalErrorLine
