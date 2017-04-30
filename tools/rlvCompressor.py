#!/usr/bin/env python2
#
# Compresses an input RLV (ReLU Verifier) instance in the sense that
# superfluous linear nodes (that are only input to one ReLU node with weight 1)
# are removed (by merging them into the ReLu node).
#
# This makes LP solving a bit easier (due to the variable number reduction) and speeds up the whole process.
#
# TODO: FLOATS need to be printed with maximal precision
import os, sys

if len(sys.argv)<2:
    print >>sys.stderr, "Error: Expected input file"
    sys.exit(1)
    
with open(sys.argv[1],"r") as inputFile:
    inputFileLines = [a.strip() for a in inputFile.readlines()]
    

# While linear nodes are there?
linearNodes = set([])
linearDefinitionLines = {}
for line in inputFileLines:
    line = line.split(" ")
    if line[0]=="Linear":
        linearNodes.add(line[1])
        linearDefinitionLines[line[1]] = (float(line[2]),line[3:])
        
# Which linear nodes are used outside of ReLUs with Weight 1?
linearNontriviallyUsedNodes = set([])
for line in inputFileLines:
    line = line.split(" ")
    if line[1].startswith("outX"):
        for i,a in enumerate(line):
            if a in linearNodes:
                    linearNontriviallyUsedNodes.add(a)
    if line[0]=="MaxPool":
        for i,a in enumerate(line):
            if (i>1) and (a in linearNodes):
                linearNontriviallyUsedNodes.add(a)
    else:
        for i,a in enumerate(line):
            if (i>1) and (a in linearNodes):
                value = float(line[i-1])
                if value!=1.0 or line[0]!="ReLU":
                    linearNontriviallyUsedNodes.add(a)
                    
# Now output the new instance, shortened by all Linear nodes that are pushed into ReLUs
for line in inputFileLines:
    line = line.split(" ")
    if line[0]=="Linear":
        if line[1] in linearNontriviallyUsedNodes:
            print " ".join(line)
        else:
            # Remove this one!
            pass
    elif line[0]=="ReLU":
        # Build the new line
        constantValue = float(line[2])
        monomials = {}
        for i in range(3,len(line),2):
            factor = float(line[i])
            variable = line[i+1]
            if (variable in linearNodes) and not (variable in linearNontriviallyUsedNodes):
                # Add the ReLU inputs
                (addedConstant,addedMonomials) = linearDefinitionLines[variable]
                constantValue += addedConstant
                for j in range(0,len(addedMonomials),2):
                    if addedMonomials[j+1] in monomials:
                        monomials[addedMonomials[j+1]] += float(addedMonomials[j])
                    else:
                        monomials[addedMonomials[j+1]] = float(addedMonomials[j])
            else:
                if variable in monomials:
                    monomials[variable] += float(factor)
                else:
                    monomials[variable] = float(factor)
        print line[0]+" "+line[1]+" "+str(constantValue)+" "+" ".join([str(b)+" "+a for (a,b) in monomials.iteritems()])
    elif line[0]=="Assert" or line[0]=="Input" or line[0]=="MaxPool":
        print " ".join(line)
    else:
        print >>sys.stderr, "Error: Unknown line type: "+line[0]
        sys.exit(1)


