#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
import string 

print 
print 'START'
print 

import argparse
parser = argparse.ArgumentParser(usage="mergeFiles.py [options]",description="Submit the jobs",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--inputdir",help='name of the directory',required=True)
parser.add_argument("--outputdir",help='name of the directory')
parser.add_argument("--dryrun",action='store_true', help='Do not submit the jobs')
args = parser.parse_args()

print "Get List of different datasets:" 
filelist = os.listdir(args.inputdir)
filelist.sort()
print filelist
datasetlist = {}
filecount = 0
for i,f in enumerate(filelist): 
    current = f.rstrip(".root").rstrip(string.digits).rstrip("_")
    if i==0: prevfile = current
    filecount+=1
    if prevfile != current or filecount==len(filelist):
        datasetlist[prevfile]=filecount
        prevfile = current
        filecount=0
        print "    "+prevfile

if len(datasetlist)==0:
    print datasetlist

outputdir=args.inputdir.rstrip('/')
if args.outputdir is not None:
    outputdir=args.outputdir
print
print "Merging inside: "+outputdir
print
for d,num in datasetlist.items():
    cmd = "hadd {idir}_{data}.root".format(idir=outputdir,data=d)
    for n in range(1,num+1):
        cmd+=" {idir}/{data}_{num}.root".format(idir=args.inputdir,data=d,num=n)
    if args.dryrun: 
        os.system("echo %s" %cmd)
    else: 
        os.system("echo %s" %cmd)
        os.system(cmd)
        
print 
print "DONE"
print


