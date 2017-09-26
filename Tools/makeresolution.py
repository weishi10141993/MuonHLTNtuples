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
YMIN = 0.01
YMAX = 1.02

ROOT.gROOT.LoadMacro("cruijff.C+")

def makePlots(inputfiles,draw,title,legends,outname):
    style = setTDRStyle()
    ROOT.gStyle.SetTitleYOffset(1.45)
    ROOT.gStyle.SetTitleXOffset(1.45)
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetStatX(.9)
    ROOT.gStyle.SetStatY(.9)

#    ROOT.TH1.AddDirectory(ROOT.kFALSE)

    #make canvas to save plots to    
    _files=[]
    for i in inputfiles:
        print "Opening %s..." %(i)
        _files.append(ROOT.TFile(i))
       
    _hists=[]
    _hist = ROOT.TH1F()
    _hist.SetDirectory(0)
    if logY: outname+="_log"    
    for i,f in enumerate(_files):
        histname = "%s" %(draw[0])
        _hist = f.Get(histname).Clone()
        _hist.SetDirectory(0)
        _hists.append(_hist)

    leg = ROOT.TLegend(0.18,0.76,0.55,0.90);
    leg.SetLineColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0)
    

    c1 = ROOT.TCanvas('c1')

    o = open(outname+'.txt', 'w')
    for i,hist in enumerate(_hists):
        if "h_" in draw[0]:
            hist.Scale(1./hist.GetEntries())
        hist.SetLineWidth(2)
        hist.SetLineColor(i+1)
        hist.SetMarkerColor(i+1)
        hist.SetTitle(title)

        YMAX=1.25*hist.GetMaximum()
        hist.GetYaxis().SetRangeUser(0.0,YMAX)

        (mean,std) = doFit(hist,outname,legends[i],i)
        o.write('%s \t %5.4f +/- %5.4f\n' %(legends[i],mean,std))
        
        c1.cd()
        if logY: c1.SetLogy()
        if logX: c1.SetLogx()

        if i==0: 
            hist.Draw("HIST")
        else:
            hist.Draw("HIST SAME")
            
        leg.AddEntry(hist,legends[i],"l")
    
    o.close()

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.04)
    latex.SetNDC(True)
    latexCMS = ROOT.TLatex()
    latexCMS.SetTextFont(61)
    latexCMS.SetTextSize(0.055)
    latexCMS.SetNDC(True)
    latexCMSExtra = ROOT.TLatex()
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextSize(0.03)
    latexCMSExtra.SetNDC(True)    
    latex.DrawLatex(0.95, 0.96, "(13 TeV)")
    
    cmsExtra = "Private" 
    latexCMS.DrawLatex(0.78,0.88,"CMS")
    yLabelPos = 0.84
    latexCMSExtra.DrawLatex(0.78,yLabelPos,"%s"%(cmsExtra))

    if len(_hists)>1:
        leg.Draw("SAME")
    
    c1.SaveAs(outname+".root");
    c1.SaveAs(outname+".pdf");
    c1.SaveAs(outname+".png");
  
def doFit(hist,outname,legend,index):
    c2 = ROOT.TCanvas("c2","c2",700,700)
    c2.cd()
    
    fit_min = hist.GetMean() - 1.25*hist.GetRMS() 
    fit_max = hist.GetMean() + 1.25*hist.GetRMS()

    gaus = ROOT.TF1("gaus","gaus",fit_min,fit_max)
    gaus.SetLineColor(4)
    gaus.SetParameters(0,hist.GetMean(),hist.GetRMS())
    hist.Fit("gaus","M0R+")

    hist.Draw("E")

    funct = ROOT.TF1()
    funct = ROOT.TF1("cruijff",ROOT.cruijff,fit_min,fit_max,5)
    funct.SetParameters(gaus.GetParameter(0), gaus.GetParameter(1), gaus.GetParameter(2), 0., 0.) #15, 0.001)             
    funct.SetParNames("Constant","Mean","Sigma", "AlphaL","AlphaR")        
    funct.SetLineColor(ROOT.kBlue)
    funct.SetLineWidth(2)
    hist.Fit("cruijff","M0R+")
    funct.Draw("SAME")

#    leg = ROOT.TLegend(0.55,0.65,0.90,0.80);
#    leg.SetLineColor(0);
#    leg.SetFillStyle(0);
#    leg.SetBorderSize(0)
#    leg.AddEntry(hist,legend,"l")
#    leg.Draw("SAME")

    latex = ROOT.TLatex()
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.04)
    latex.SetNDC(True)
    latexCMS = ROOT.TLatex()
    latexCMS.SetTextFont(61)
    latexCMS.SetTextSize(0.055)
    latexCMS.SetNDC(True)
    latexCMSExtra = ROOT.TLatex()
    latexCMSExtra.SetTextFont(52)
    latexCMSExtra.SetTextSize(0.03)
    latexCMSExtra.SetNDC(True)    
    latex.DrawLatex(0.95, 0.96, "(13 TeV)")
    
    cmsExtra = "Private" 
    latexCMS.DrawLatex(0.78,0.88,"CMS")
    yLabelPos = 0.84
    latexCMSExtra.DrawLatex(0.78,yLabelPos,"%s"%(cmsExtra))

    latexFit = ROOT.TLatex()
    latexFit.SetTextFont(61)
    latexFit.SetTextSize(0.035)
    latexFit.SetNDC(True)
    latexFit.DrawLatex(0.19,0.88, legend)
    
    latexFit1 = ROOT.TLatex()
    latexFit1.SetTextFont(42)
    latexFit1.SetTextSize(0.035)
    latexFit1.SetNDC(True)
    for par in range(funct.GetNpar()-3):
        yPos = 0.84-0.04*(float(par))
        latexFit1.DrawLatex(0.19, yPos, "%s = %5.2g #pm %5.2g"%(funct.GetParName(par+1),funct.GetParameter(par+1),funct.GetParError(par+1))) 
    
    saveas = "%s_h%s_fit" %(outname,index)
    c2.SaveAs(saveas+".root");
    c2.SaveAs(saveas+".pdf");
    c2.SaveAs(saveas+".png");

    mean = gaus.GetParameter(1)
    sig  = gaus.GetParameter(2)

    return mean,sig


#### ========= MAIN =======================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="makeEffs.py [options]",description="Extract Efficiencies from file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i","--inputfiles", type=str, help='The list of input files, comma separated if more than one file',required=True,nargs=1)
    parser.add_argument("--draw",type=str, help='the plot to be compared',nargs=1)
    parser.add_argument("--title",dest="title", type=str, help='title of the histo;xaxis;yaxis')    
    parser.add_argument("--leg", type=str,help='the list of legends to be used, comma separated',required=True,nargs=1)
    parser.add_argument("--outdir",dest="fdir", default="test/", help='name of the outputfile')
    parser.add_argument("--outname",dest="outname", default="efficiency", help='name of the outputfile')
    parser.add_argument("--logy", action='store_true', help='activate LogY')
    parser.add_argument("--logx", action='store_true', help='activate LogX')
    parser.add_argument("--yrange", type=str,help='lower and upper y limit', nargs=1)

#    parser.add_argument("-o","--ofolder",dest="output", default="plots/", help='folder name tose store results')
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

    if args.draw is not None:
        draw = args.draw[0].split(",")
        makePlots(files,draw,args.title,legends,args.fdir+args.outname)

