#!/usr/bin/env python2
# Check how often learning leads to 100% accuracy
import os, sys, subprocess

for i in range(0,100):
    experiment = subprocess.Popen("./train.sh",shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    acc = None
    for line in experiment.stdout.readlines():
        if line.find("Test net output #0: accuracy =")!=-1:
            acc = line.split(" ")[-1]
    with open("poolingAccurracies.txt","a") as outFile:
        outFile.write(acc)

