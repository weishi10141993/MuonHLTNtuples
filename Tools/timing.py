# import ROOT in batch mode
import sys,getopt,os

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
import math
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from setTDRStyle import setTDRStyle


def getTiming(inputfile,pathname,beginmodule,endmodule):
    
    print "Opening %s..." %(inputfile)
    _file=ROOT.TFile(inputfile)

    
    basename = "DQMData/Run 999999/HLT/Run summary/TimerService/process TEST paths/path %s" %(pathname)
    histname = basename+"/module_time_real_total"
    timename = basename+"/path time_real"
    print histname
    _time = _file.Get(timename).Clone()
    totaltime = _time.GetMean()
    totalnum  = _time.GetEntries()
    _hist = _file.Get(histname).Clone()
    
    xlow = 0.
    for bin in xrange(_hist.GetNbinsX()):
        if beginmodule == _hist.GetXaxis().GetBinLabel(bin): 
            xlow = bin
        if endmodule == _hist.GetXaxis().GetBinLabel(bin):
            xup = bin
    
    print
    print "Timing: " 
    for bin in xrange(xlow,xup+1):
        print "   {module:55s} = {time:5.3f} ms".format(module=_hist.GetXaxis().GetBinLabel(bin),time=_hist.GetBinContent(bin)/totalnum)


    print
    print "Total Timing: {time:5.3f}".format(time=totaltime) 
    print "   {bmodule} - {emodule} = {time:5.3f} ms".format(bmodule=beginmodule,emodule=endmodule,time=_hist.Integral(xlow,xup)/totalnum)
    print
    
def printModules(inputfile,pathname):
    print "Opening %s..." %(inputfile)
    _file=ROOT.TFile(inputfile)

    basename = "DQMData/Run 999999/HLT/Run summary/TimerService/process TEST paths/path %s" %(pathname)
    histname = basename+"/module_time_real_total"
    print histname
    _hist = _file.Get(histname).Clone()
    
    print
    print "Printing modules..."
    for bin in xrange(_hist.GetNbinsX()):
        print "   {module:55s}".format(module=_hist.GetXaxis().GetBinLabel(bin))
    print
                                       


#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="makeEffs.py [options]",description="Extract Efficiencies from file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i","--inputfile", type=str, help='The list of input files, comma separated if more than one file',required=True)
    parser.add_argument("--path",type=str, help='Path name',default="HLT_TestMuIterL3_Mu50_v1")
    parser.add_argument("--time",type=str, help='the numerator,denominator to be used', nargs=1)
    parser.add_argument("--printmodules", action='store_true', help='printAllModules')

#    parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name tose store results')
    args = parser.parse_args()
    
    if args.printmodules:
        printModules(args.inputfile,args.path)
        sys.exit()

    if args.time is not None: 
        modules = args.time[0].split(",")
        getTiming(args.inputfile,args.path,modules[0],modules[1])
