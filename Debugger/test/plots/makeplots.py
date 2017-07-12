#!/usr/bin/env python
import sys,os
import re

ODIR=sys.argv[2]  # "/afs/cern.ch/user/f/folguera/www/private/MUO/170401_MTV/ZMM/iter02/"
IFILE=sys.argv[1] #"170403_debuggingIO/DQM_90X_iter02_L2ROI_TrajParams_ZMM.root"
#DOWHAT="DQMCombo"
#DOWHAT="DQM"
#DOWHAT="DeltaRCombo"
DOWHAT=sys.argv[3]

tracklistB=[
    "hltIterL3MuonSeededOutIn",                     ##IterOI 
    #IO            "hltIterL3MuonPixel",                           ##IterIO
    #IO            "hltIter0IterL3MuonCtfWithMaterial",            ##IterIO
    #IO            "hltIter0IterL3MuonTrackSelectionHighPurity",   ##IterIO
    #IO            "hltIter2IterL3MuonCtfWithMaterial",            ##IterIO
    #IO            "hltIter2IterL3MuonTrackSelectionHighPurity",   ##IterIO
    #IO            "hltIter2IterL3MuonMerged",                     ##IterIO 
    "hltIterL3MuonMerged",
            ]  

tracklistA=[
    "hltL3TkFromL2OIState",                     ##Cascade
#   "hltL3TkFromL2OIHit",                       ##Cascade
#   "hltL3TrackCandidateFromL2",                ##Cascade
    #"hltIter0HighPtTkMuPixel",                      ##TkMu 
    #"hltIter0HighPtTkMuCtfWithMaterial",            ##TkMu
    #"hltIter0HighPtTkMuTrackSelectionHighPurity",   ##TkMu
    #"hltIter2HighPtTkMuCtfWithMaterial",            ##TkMu
    #"hltIter2HighPtTkMuTrackSelectionHighPurity",   ##TkMu
    #"hltIter2HighPtTkMuMerged",                     ##TkMu
    ]

tracklistCombo2016=[
    "hltL3TkMergeStep1",
    "hltIter2HighPtTkMuMerged",  
    "hltTkFromL2AndTkMu"
    ]
tracklistCombo2017=[
    "hltIterL3OIMuCtfWithMaterial",
    "hltIterL3MuonMerged",
    "hltIterL3MuonAndMuonFromL1Merged",
    ]

plotsDQM=["effic",
          "efficPt",
          "effic_vs_dz",
          "effic_vs_dxy",
          "effic_vs_hit",
          "effic_vs_layer",
          "fakerate",
          "fakeratePt",
          "fakerate_vs_dz",
          "fakerate_vs_hit",
          "fakerate_vs_layer",
          "duplicatesRate",
          "duplicatesRate_Pt",
          "duplicatesRate_chi2",
          ]

plotsEff=["effL3L2Eta",
          "effL3L2Pt",
          "effL3L2Phi",
          "effL3L2NPV",
          ]

if __name__ == '__main__':
    
    if "DQM" in DOWHAT: 
        for plot in plotsDQM:
            if "Combo" in DOWHAT: 
                if "Displaced" in IFILE:
                    yrange="--yrange='0.01,1.0'"
                else:
                    yrange="--yrange='0.80,1.0'"
                
                if 'duplicate' in plot: yrange="--yrange='0.01,0.1'"
                if 'fake' in plot: yrange="--yrange='0.01,0.1'"
                if 'Pt' in plot: yrange+=" --logx"
                
                if len(IFILE.split(','))==2:
                    for i,t in enumerate(tracklistCombo2017):
                        legA="2016"
                        legB="2017"
                        plotA="{fld}_hltAssociatorByHits/{plot}".format(fld=tracklistCombo2016[i],plot=plot)
                        plotB="{fld}_hltAssociatorByHits/{plot}".format(fld=tracklistCombo2017[i],plot=plot)
                        name="{trk}_{plot}".format(trk=tracklistCombo2017[i].replace("hlt",""),plot=plot)
                        cmd = "python plotter.py --outdir {odir} --inputfile {ifile} --outname {name} --leg '{legA},{legB}' --dqm '{plotA},{plotB}' {yrange}".format(odir = ODIR, ifile=IFILE, name=name,legA=legA,legB=legB,plotA=plotA,plotB=plotB,yrange=yrange )
                        print "Running... "+ cmd
                        os.system(cmd)
                else: 
                    legA="TkMu (2016)"
                    legB="Cascade OR TkMu (2016)"
                    for t in tracklistCombo2017:
                        legC="" 
                        if "hltIterL3OIMuCtfWithMaterial" in t: continue
                        if "hltIterL3MuonMerged" in t: 
                            legC="IterL3: OI+IOL2"
                        if "hltIterL3MuonAndMuonFromL1Merged" in t: 
                            legC="IterL3: OI+IOL2+IOL1"                       
                            
                        plotA="{fld}_hltAssociatorByHits/{plot}".format(fld="hltIter2HighPtTkMuMerged",plot=plot)
                        plotB="{fld}_hltAssociatorByHits/{plot}".format(fld="hltTkFromL2AndTkMu",plot=plot)
                        plotC="{fld}_hltAssociatorByHits/{plot}".format(fld=t,plot=plot)
                        name="{trk}_Combo_{plot}".format(trk=t,plot=plot)
                        cmd = "python plotter.py --outdir {odir} --inputfile {ifile} --outname {name} --leg '{legA},{legB},{legC}' --dqm '{plotA},{plotB},{plotC}' {yrange}".format(odir = ODIR, ifile=IFILE, name=name,legA=legA,legB=legB,legC=legC,plotA=plotA,plotB=plotB,plotC=plotC,yrange=yrange )
                        print "Running... "+ cmd
                        os.system(cmd)
                
            else: 
                for i,t in enumerate(tracklistA):
                    yrange="--yrange='0.0,1.0'"
                    if 'duplicate' in plot: yrange="--yrange='0.0,0.1'"
                    if 'fake' in plot: yrange="--yrange='0.0,0.1'"
                    if 'Pt' in plot: yrange="--logx"

                    legA="2016"
                    legB="2017"
                    plotA="{fld}_hltAssociatorByHits/{plot}".format(fld=tracklistA[i],plot=plot)
                    plotB="{fld}_hltAssociatorByHits/{plot}".format(fld=tracklistB[i],plot=plot)
                    name="{trk}_{plot}".format(trk=tracklistA[i].replace("HighPtTkMu",""),plot=plot)
                    cmd = "python plotter.py --outdir {odir} --inputfile {ifile} --outname {name} --leg '{legA},{legB}' --dqm '{plotA},{plotB}' {yrange}".format(odir = ODIR, ifile=IFILE, name=name,legA=legA,legB=legB,plotA=plotA,plotB=plotB,yrange=yrange )
                    print "Running... "+ cmd
                    os.system(cmd)

            
    if "DeltaR" in DOWHAT:
        for plot in plotsEff:
            if "Combo" in DOWHAT:
                if "Displaced" in IFILE:
                    yrange="--yrange='0.01,1.0'"
                else:
                    yrange="--yrange='0.85,1.0'"
                legA="IterL3 (OI)"
                legB="IterL3 (IO-L2)"
                legC="IterL3 (Combo)"
                legD="IterL3 (Combo w IO-L1)"
                name="{plot}".format(plot=plot)
                cmd = "python plotter.py --outdir {odir} --inputfile '{ifile}' --outname {name} --leg '{legA},{legB},{legC},{legD}' --title 'L3/L2 efficiency' --eff '{plot}_num,{plot}_den' {yrange}".format(odir = ODIR, ifile=IFILE, name=name,legA=legA,legB=legB,legC=legC,legD=legD,plot=plot,yrange=yrange )
                print "Running... "+ cmd
                os.system(cmd)

##                yrange="--yrange='0.85,1.0'"
##                legA="Combo (2017)"
##                legB="OI (2017)"
##                legC="IO (2017)"
##                name="{plot}".format(plot=plot)
##                cmd = "python plotter.py --outdir {odir} --inputfile '{ifile}' --outname {name} --leg '{legA},{legB},{legC}' --title 'L3/L2 efficiency' --eff '{plot}_num,{plot}_den' {yrange}".format(odir = ODIR, ifile=IFILE, name=name,legA=legA,legB=legB,legC=legC,plot=plot,yrange=yrange )
            else:
                if "Displaced" in IFILE:
                    yrange="--yrange='0.01,1.0'"
                else:
                    yrange="--yrange='0.85,1.0'"

                legA="2016"
                legB="2017"
                name="{plot}".format(plot=plot)
                cmd = "python plotter.py --outdir {odir} --inputfile '{ifile}' --outname {name} --leg '{legA},{legB}' --title 'L3/L2 efficiency' --eff '{plot}_num,{plot}_den' {yrange}".format(odir = ODIR, ifile=IFILE, name=name,legA=legA,legB=legB,plot=plot,yrange=yrange )
                print "Running... "+ cmd
                os.system(cmd)
            
