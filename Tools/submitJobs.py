#!/usr/bin/env pyothon
import os, re
import commands
import math, time
import sys

print 
print 'START'
print 

import argparse
parser = argparse.ArgumentParser(usage="submitJobs.py [options]",description="Submit the jobs",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--inputFiles", type=str, help='The list of input files')
parser.add_argument("--inputSecondaryFiles", type=str, help='The list of Secondary input files')
parser.add_argument("--useSecondaryFiles", action='store_true', help='Use secondary Files')
parser.add_argument("--cfg",type=str, help='The config file')
parser.add_argument("--outputFile", help='prefix of the outputfile')
parser.add_argument("--secondaryOutput", help='prefix of the Secondary outputfile')
parser.add_argument("--outdir", default="test/", help='name of the outputDir')
parser.add_argument("--queue", default="workday", help='queue name:  espresso (20 min), microcentury (1h), longlunch (2h), workday (8h), tomorrow (1d)')
parser.add_argument("--clean",action='store_true', help='Clean all the working folders')
parser.add_argument("--dryrun",action='store_true', help='Do not submit the jobs')
parser.add_argument("--nfilesperjob",type=int, help="Number of files per job, if nothing is specified one file per job will be run")

args = parser.parse_args()

if args.clean:
    print "Cleaning..."
    print "  - Cleaning executables"
    os.system("rm -rf exec")
    print "  - Cleaning logs"
    os.system("rm -rf batchlogs")
    print 
    print "Done!"
    print 
    sys.exit()
else: 
    if args.inputFiles is None: 
        print "You need to provide a --inputFiles, I don't know what to do otherwise"
        print 
        sys.exit()
    if args.cfg is None:
        print "You need to provide a --cfg, I don't know what to do otherwise"
        print 
        sys.exit()


########   customization  area #########
OutputDir = args.outdir
## outputfile name is determined from the List of files: 
#for i in $(cat AllDatasets_RECO.txt ); 
#   do output="${i//\//_}"; 
#   echo "saving in List$output.txt"; 
#   das_client --query="file dataset=${i}" --limit=1000 | grep .root &> List${output}.txt &  
#done
## Example ListName: List_RelValZMM_13_CMSSW_9_2_0-PU25ns_91X_upgrade2017_realistic_v5-v1_GEN-SIM-RECO.txt
OutputFileNames = "" 
SecondaryOutputFileNames = ""
ScriptNames = ""
if args.outputFile is None:
    target = args.inputFiles
    print target
    rdict = {
        'fileList/List_': '',
        'fileList/': '',
        '_CMSSW_9_2_0-PU25ns_91X_upgrade2017_realistic_v5': '',
        '_CMSSW_9_2_0-91X_upgrade2017_realistic_v5': '',
        '-v1_GEN-SIM-RECO': '',
        '-v1_GEN-SIM-DIGI-RAW': '',
        '-v1_RAW': '',
        '-v1_AOD': '',
        '_RAW-RECO': '',
        '.txt': ''
        }
    robj = re.compile('|'.join(rdict.keys()))
    outname = robj.sub(lambda m: rdict[m.group(0)], target)
    print outname
    OutputFileNames = OutputDir+"muonNtuple_"+outname # base of the output file name, they will be saved in res directory
    SecondaryOutputFileNames = OutputDir+"TP_"+outname # base of the output file name, they will be saved in res directory
    ScriptNames = outname
else:
    OutputFileNames = OutputDir+args.outputFile
    SecondaryOutputFileNames = OutputDir+args.secondaryOutput
    ScriptNames = args.outputFile



CfgName = args.cfg # script to be used with cmsRun
ScriptFolder = CfgName.replace(".py","")
FileListRECO = args.inputFiles # list with all the file directories

print "Getting list of secondary files..." 
FileListRAW = "" # list with all the file directories
if args.useSecondaryFiles and args.inputSecondaryFiles is None:
    FileListRAW = args.inputFiles.replace("GEN-SIM-RECO","GEN-SIM-DIGI-RAW")
    FileListRAW = FileListRAW.replace("AOD","RAW")
if args.inputSecondaryFiles is not None: 
    FileListRAW  =  args.inputSecondaryFiles   
print FileListRAW

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

with open(FileListRECO) as f:
    for i, l in enumerate(f):
        pass
    NumberOfFiles=i+1

NumberOfFilesPerJob=1
if args.nfilesperjob is not None:
    NumberOfFilesPerJob=args.nfilesperjob

queue = args.queue  # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 
########   customization end   #########

path = os.getcwd()
print 'do not worry about folder creation:'
if not os.path.exists("exec"):  os.makedirs("exec")
if not os.path.exists("exec/%s" %(ScriptFolder)):  os.makedirs("exec/%s" %(ScriptFolder))
if not os.path.exists("exec/%s/%s" %(ScriptFolder,ScriptNames)):  os.makedirs("exec/%s/%s" %(ScriptFolder,ScriptNames))
if not os.path.exists("batchlogs"):  os.makedirs("batchlogs")
if not os.path.exists(OutputDir):    os.makedirs(OutputDir)
print 

##### loop for creating and sending jobs #####
print "Creating jobs for %s" %FileListRECO
import itertools
NumberOfJobs=0
with open(FileListRECO) as T:
    for i, sli in enumerate(iter(lambda:list(itertools.islice(T, NumberOfFilesPerJob)), []), 0):
        with open("exec/{fld}/{subfld}/list_{n}.txt".format(fld=ScriptFolder,subfld=ScriptNames,n=i), "w") as f:
            f.writelines(sli)
            NumberOfJobs+=1

print "Number of jobs created: %s" %NumberOfJobs

for x in range(0, int(NumberOfJobs)):
    ##    with open(FileListRECO) as f:
    ##### creates directory and file list for job #######
                
#    os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' "+FileListRECO+" > exec/%s/%s/list_%s.txt" %(ScriptFolder,ScriptNames,x))
#    inputfiles = "exec/%s/%s/list_%s.txt" %(ScriptFolder,ScriptNames,x)
#    scriptname = 'exec/%s/%s/job_%s.sh' %(ScriptFolder,ScriptNames,x)
#    os.system("sed '"+str(1+interval*(x-1))+","+str(interval*x)+"!d' "+FileListRECO+" > exec/%s/%s/list_%s.txt" %(ScriptFolder,ScriptNames,x-1))
    inputfiles = "exec/%s/%s/list_%s.txt" %(ScriptFolder,ScriptNames,x)
    scriptname = 'exec/%s/%s/job_%s.sh' %(ScriptFolder,ScriptNames,x)

    print " ---> %s " %(scriptname)
    ##### creates jobs #######
    with open(scriptname, 'w') as fout:
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("echo 'WORKDIR ' ${PWD}\n")
        fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
        fout.write("cd "+str(path)+"\n")
        fout.write("cmsenv\n")
        ### defining variables: 
        outfilename = OutputFileNames+"_"+str(x-1)+".root"
        secfilename = SecondaryOutputFileNames+"_"+str(x-1)+".root"
        if args.useSecondaryFiles:
            fout.write("cmsRun "+CfgName+" outputFile='"+outfilename+"' secondaryOutputFile='"+secfilename+"' inputFiles_clear inputFiles_load='"+inputfiles+"' secondaryInputFiles_clear secondaryInputFiles_load='"+FileListRAW+"'\n")                
        else: 
            fout.write("cmsRun "+CfgName+" outputFile='"+outfilename+"' secondaryOutputFile='"+secfilename+"' inputFiles_clear inputFiles_load='"+inputfiles+"'\n")                
        fout.write("echo 'STOP---------------'\n")
        fout.write("echo\n")
        fout.write("echo\n")
        os.system("chmod 755 %s" %scriptname)
            
###### create submit.sub file ####
with open('submit.sub', 'w') as fout:
    fout.write("executable              = $(filename)\n")
    fout.write("arguments               = $(ClusterId)$(ProcId)\n")
    fout.write("output                  = batchlogs/%s_%s$(ClusterId).$(ProcId).out\n" %(ScriptFolder,ScriptNames))
    fout.write("error                   = batchlogs/%s_%s$(ClusterId).$(ProcId).err\n" %(ScriptFolder,ScriptNames))
    fout.write("log                     = batchlogs/%s_%s$(ClusterId).log\n"%(ScriptFolder,ScriptNames))
    fout.write('+JobFlavour = "%s"\n' %(queue))
    fout.write("\n")
    fout.write("queue filename matching (exec/%s/%s/job_*.sh)\n" %(ScriptFolder,ScriptNames))
    
###### sends bjobs ######
os.system("cat submit.sub")
if not args.dryrun: 
    os.system("condor_submit submit.sub")
   
print
print "your jobs:"
os.system("condor_q")
print
print 'END'
print
