# import ROOT in batch mode
import sys,getopt,os

oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
import ROOT
import math
ROOT.gROOT.SetBatch(True)
sys.argv = oldargv

from setTDRStyle import setTDRStyle

logY = False
logX = False
YMIN = 0.5
YMAX = 1.02

def makePlots(inputfiles,draw,title,legends,outname):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.Reset()
    #make canvas to save plots to
    c1 = ROOT.TCanvas('c1')
    
    _files=[]
    for i in inputfiles:
        _files.append(ROOT.TFile(i))
    
    _hists=[]
    _hist = ROOT.TH1F()
    _hist.SetDirectory(0)
    ROOT.TH1.AddDirectory(ROOT.kFALSE)   

    if logY: outname+="_log"
    o = open(outname+'.txt', 'w')
    for i,f in enumerate(_files):
        for k,p in enumerate(draw):
            histname = "muonDebugger/%s" %(p)
            _hist = f.Get(histname).Clone()
            _hists.append(_hist)
            if (len(draw)>1): 
                o.write('%s \t %5.4f \t %5.4f\n' %(legends[k],_hist.GetMean(),_hist.GetRMS()))
            else:
                o.write('%s \t %5.4f \t %5.4f\n' %(legends[i],_hist.GetMean(),_hist.GetRMS()))
    o.close()

    leg = ROOT.TLegend(0.70,0.65,0.90,0.80);
    leg.SetLineColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0)
    
    k=0
    for hist in _hists:
        hist.Scale(1/hist.GetEntries())
        hist.SetLineWidth(2)
        hist.SetLineColor(k+1)
        hist.SetTitle(title)
        if k==0 :
            if logY: c1.SetLogy()
            if logX: c1.SetLogx()
            if "Vs" in draw[0]:
                hist.Draw("BOX")
            else: 
                hist.Draw()
        else:
            hist.Draw("SAME")

        leg.AddEntry(hist,legends[k],"l")
        k+=1
    
    if len(_hists)>1:
        leg.Draw("SAME")
    
    c1.SaveAs(outname+".root");
    c1.SaveAs(outname+".pdf");
    c1.SaveAs(outname+".png");


def makePlotsFromDQM(inputfiles,draw,legends,outname):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.Reset()
    c1 = ROOT.TCanvas('c1')
    
    _files=[]
    for i in inputfiles:
        _files.append(ROOT.TFile(i))
        
    _hists=[]
    _hist = ROOT.TH1F()
    _hist.SetDirectory(0)
#    ROOT.TH1.AddDirectory(ROOT.kFALSE)   

    if logY: outname+="_log"
    o = open(outname+'.txt', 'w')
    if len(_files)!=len(draw):
        for i,f in enumerate(_files):
            for k,p in enumerate(draw):
                histname = "DQMData/Run 1/HLT/Run summary/Tracking/ValidationWRTtp/%s" %(p)
                _hist = f.Get(histname).Clone()
                
                ## get mean
                mean = 0
                err = 0 
                for b in xrange(_hist.GetNbinsX()): 
                    mean = mean+_hist.GetBinContent(b)
                    err  = err+_hist.GetBinError(b)*_hist.GetBinError(b)
 
                mean = mean/_hist.GetNbinsX() 
                err  = math.sqrt(err/_hist.GetNbinsX())

                _hists.append(_hist)
                if (len(draw)>1): 
                    o.write('%s \t %5.3f +/- %5.3f\n' %(legends[k],mean,err))
                else:
                    o.write('%s \t %5.3f +/- %5.3f\n' %(legends[i],mean,err))
    else: 
        for i,f in enumerate(_files):
            histname = "DQMData/Run 1/HLT/Run summary/Tracking/ValidationWRTtp/%s" %(draw[i])
            print histname
            _hist = f.Get(histname).Clone()
            
            ## get mean
            mean = 0
            err = 0
            for b in xrange(_hist.GetNbinsX()): 
                mean = mean+_hist.GetBinContent(b)
                err  = err+_hist.GetBinError(b)*_hist.GetBinError(b)
            mean = mean/_hist.GetNbinsX() 
            err  = math.sqrt(err/_hist.GetNbinsX())
            
            _hists.append(_hist)
            o.write('%s \t %5.3f +/- %5.3f\n' %(legends[i],mean,err))
            
            
    o.close()
    leg = ROOT.TLegend(0.30,0.15,0.70,0.3);
    if "Displaced" in outname:  leg = ROOT.TLegend(0.60,0.70,0.90,0.85);
    leg.SetLineColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0)
    
    k=0
    for hist in _hists:
        hist.SetLineWidth(2)
        hist.SetLineColor(k+1)
#        hist.SetTitle(title)
        hist.GetYaxis().SetRangeUser(YMIN,YMAX)

        if k==0 :
            if logY: c1.SetLogy()
            if logX: c1.SetLogx()
            if "Vs" in draw[0]:
                hist.Draw("BOX")
            else: 
                hist.Draw()
        else:
            hist.Draw("SAME")

        leg.AddEntry(hist,legends[k],"l")
        k+=1
    
    if len(_hists)>1:
        leg.Draw("SAME")
    
    c1.SaveAs(outname+".root");
    c1.SaveAs(outname+".pdf");
    c1.SaveAs(outname+".png");



def makeEfficiencies(inputfiles,num,den,title,legends,outname):
    ROOT.gStyle.SetOptStat(0)
    
    ROOT.gROOT.Reset()
    #make canvas to save plots to
    c = ROOT.TCanvas('c')

    _files=[]
    for i in inputfiles:
        _files.append(ROOT.TFile(i))
    
    _effs=[]
    _num = ROOT.TH1F()
    _den = ROOT.TH1F()
    _num.SetDirectory(0)
    _den.SetDirectory(0)
#    ROOT.TH1.AddDirectory(ROOT.kFALSE)   

    o = open(outname+'.txt', 'w')
    for i,f in enumerate(_files):
        numname = "muonDebugger/%s" %(num)
        denname = "muonDebugger/%s" %(den)
        print "Trying to get %s and %s from %s" %(numname,denname, f.GetName())
        _num = f.Get(numname).Clone()
        _den = f.Get(denname).Clone()
        o.write('%s \t %5.4f\n' %(legends[i],_num.GetEntries()/_den.GetEntries()))
 
        eff = ROOT.TEfficiency(_num,_den)
        _effs.append(eff)
    o.close()

    leg = ROOT.TLegend(0.30,0.15,0.70,0.3);
#    if "dxy" in outname:  leg = ROOT.TLegend(0.60,0.70,0.90,0.85);
#    if "dz" in outname:  leg = ROOT.TLegend(0.60,0.70,0.90,0.85);
#    if "Chi2" in outname:  leg = ROOT.TLegend(0.60,0.70,0.90,0.85);
    if "Displaced" in outname:  leg = ROOT.TLegend(0.60,0.70,0.90,0.85);
#    if "JPsi" in outname:  leg = ROOT.TLegend(0.30,0.70,0.70,0.85);
    leg.SetLineColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0)
    
    k=0
    for eff in _effs:
        if k==0 :
            eff.SetLineWidth(2)
            eff.SetLineColor(k+1)
            eff.SetTitle(title)
            eff.Draw()
            
            ## SET CORRECT AXIS VALUES
            ROOT.gPad.Update()
            graph = eff.GetPaintedGraph()
            graph.SetMinimum(YMIN)
            graph.SetMaximum(YMAX)
            if ("NPU" in outname) or ("NPV" in outname): graph.GetXaxis().SetRangeUser(20,70)
##            if ("dxy" or "dz") in outname: graph.SetMinimum(0.5)
##            if "eff2" in num: graph.SetMinimum(0.4)
##            if "JPsi" in outname: graph.SetMinimum(0.0)

            ROOT.gPad.Update()
            
        else:
            eff.SetLineWidth(2)
            eff.SetLineColor(k+1)
            eff.SetTitle(title)
            eff.Draw("SAME")
        
        leg.AddEntry(eff,legends[k],"l")
        k+=1
    
    leg.Draw("SAME")
    
    c.SaveAs(outname+".root");
    c.SaveAs(outname+".pdf");
    c.SaveAs(outname+".png");
  
    
#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="makeEffs.py [options]",description="Extract Efficiencies from file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i","--inputfiles", type=str, help='The list of input files, comma separated if more than one file',required=True,nargs=1)
    parser.add_argument("--eff",type=str, help='the numerator,denominator to be used', nargs=1)
    parser.add_argument("--draw",type=str, help='the plot to be drawn',nargs=1)
    parser.add_argument("--dqm",type=str, help='the plot to be drawn',nargs=1)
    parser.add_argument("--title",dest="title", type=str, help='title of the histo;xaxis;yaxis')    
    parser.add_argument("--leg", type=str,help='the list of legends to be used, comma separated',required=True,nargs=1)
    parser.add_argument("--outdir",dest="fdir", default="test/", help='name of the outputfile')
    parser.add_argument("--outname",dest="outname", default="efficiency", help='name of the outputfile')
    parser.add_argument("--logy", action='store_true', help='activate LogY')
    parser.add_argument("--logx", action='store_true', help='activate LogX')
    parser.add_argument("--yrange", type=str,help='lower and upper y limit', nargs=1)

#    parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name to store results')
    args = parser.parse_args()
    files = args.inputfiles[0].split(",")
    legends = args.leg[0].split(",")
    logY = args.logy
    logX = args.logx
    if args.yrange is not None:
        yranges = args.yrange[0].split(",")
        YMIN = float(yranges[0])
        YMAX = float(yranges[1])

    if not os.path.exists(args.fdir): 
        os.makedirs(args.fdir); 
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+args.fdir)
    
    if args.eff is not None: 
        effs = args.eff[0].split(",")
        makeEfficiencies(files,effs[0],effs[1],args.title,legends,args.fdir+args.outname)
    
    if args.draw is not None:
        draw = args.draw[0].split(",")
        makePlots(files,draw,args.title,legends,args.fdir+args.outname)

    if args.dqm is not None:
        dqm = args.dqm[0].split(",")
        makePlotsFromDQM(files,dqm,legends,args.fdir+args.outname)



