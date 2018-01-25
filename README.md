# MuonHLTNtuples

    cmsrel CMSSW_9_2_15  (-> update to the appropriate CMSSW version!)   
    cd CMSSW_9_2_15/src
    cmsenv    
    git cms-addpkg HLTrigger/Configuration    
    git clone git@github.com:weishi10141993/MuonHLTNtuples.git    
    cd MuonHLTNtuples/  
    git checkout -b YOUR_BRANCH_NAME    
    cd ..  
    scramv1 b   

## Produce Ntuples: 
First thing to do is to download a configuration file with the desired HLT paths, extract an HLT configuration from ConfDB using the [hltGetConfiguration](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHltGetConfiguration) script, here's an example: 

`hltGetConfiguration /online/collisions/2017/2e34/v4.2/HLT/V6 --input file:/eos/cms/store/user/folguera/ROOT/BA7BEDED-948E-E711-AFC7-02163E011FAF.root  --paths HLTriggerFirstPath,HLT_IsoMu27_v14,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v13,HLT_Mu50_v12,HLT_DoubleMu43NoFiltersNoVtx_v3,HLT_Dimuon25_Jpsi_v13,DST_DoubleMu3_noVtx_CaloScouting_v5,HLTriggerFinalPath --output none --full --offline --data --unprescale --process TEST --globaltag auto:run2_hlt_GRun > hltConfig_Ntuple.py`

Then at the end of it you should add different lines to run the Ntuple-producer: 

    process.muonNtuples = cms.EDAnalyzer("MuonNtuples",
                       offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
                       offlineMuons             = cms.InputTag("muons"),
                       triggerResult            = cms.untracked.InputTag("TriggerResults::TEST"),
                       triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::TEST"),
                       tagTriggerResult         = cms.untracked.InputTag("TriggerResults::HLT"),
                       tagTriggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                       L3Candidates             = cms.untracked.InputTag("hltIterL3MuonCandidates"),
                       L2Candidates             = cms.untracked.InputTag("hltL2MuonCandidates"),
                       L1Candidates             = cms.untracked.InputTag('hltGtStage2Digis','Muon'), 
                       TkMuCandidates           = cms.untracked.InputTag("hltIterL3OIL3MuonCandidates"),
                       L3OIMuCandidates         = cms.untracked.InputTag("hltIterL3OIL3MuonCandidates"),
                       L3IOMuCandidates         = cms.untracked.InputTag("hltIterL3IOFromL2MuonCandidates"),         
                       theTrackOI               = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity"), 
                       theTrackIOL2             = cms.untracked.InputTag("hltIter2IterL3MuonMerged"),
                       theTrackIOL1             = cms.untracked.InputTag("hltIter2IterL3FromL1MuonMerged"), 
                       lumiScalerTag            = cms.untracked.InputTag("scalersRawToDigi"),
                       puInfoTag                = cms.untracked.InputTag("addPileupInfo"),
                       genParticlesTag          = cms.untracked.InputTag("genParticles"),
                       doOffline                = cms.untracked.bool(True)
    )
    process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
    )
    process.HLTValidation = cms.EndPath(
        process.muonNtuples
    )

Then you can run the configuration file with cmsRun locally, using the batch (Tools/submitJobs.py) or using crab (i.e. `python Tools/crab_HLTMC_eff_multicrab_cfg.py`). 

## Run debugger (for event-by-event comparison) 
First thing to do is to download a configuration file with the desired HLT paths, extract an HLT configuration from ConfDB using the hltGetConfiguration script, here's an example: 

`hltGetConfiguration /online/collisions/2017/2e34/v4.2/HLT/V6 --input /store/mc/PhaseISpring17DR/SingleMu_Pt1To1000_FlatRandomOneOverPt/AODSIM/NoPUNZS_90X_upgrade2017_realistic_v20-v1/100000/0404BCB7-0234-E711-8A42-001E67DDC0FB.root  --paths HLTriggerFirstPath,HLT_IsoMu27_v14,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v13,HLT_Mu50_v12,HLT_DoubleMu43NoFiltersNoVtx_v3,HLT_Dimuon25_Jpsi_v13,DST_DoubleMu3_noVtx_CaloScouting_v5,HLTriggerFinalPath --output none --full --offline --mc --unprescale --process TEST --globaltag 92X_upgrade2017_TSG_For90XSamples_V2 --l1-emulator FullMC --l1 L1Menu_Collisions2017_v4 > hltConfig_debugger.py`

For MC, be sure for each input file above, put its parent files as secondary files in the hlt config. In this case of [AODSIM](https://cmsweb.cern.ch/das/request?input=file%3D%2Fstore%2Fmc%2FPhaseISpring17DR%2FSingleMu_Pt1To1000_FlatRandomOneOverPt%2FAODSIM%2FNoPUNZS_90X_upgrade2017_realistic_v20-v1%2F100000%2F0404BCB7-0234-E711-8A42-001E67DDC0FB.root&instance=prod%2Fglobal), its parent file are [GEN-SIM-RAW](https://cmsweb.cern.ch/das/request?input=parent%20file%3D/store/mc/PhaseISpring17DR/SingleMu_Pt1To1000_FlatRandomOneOverPt/AODSIM/NoPUNZS_90X_upgrade2017_realistic_v20-v1/100000/0404BCB7-0234-E711-8A42-001E67DDC0FB.root&instance=prod/global&idx=0&limit=10).

Then at the end of it you should add different lines to run the Debugger: 

You need to add the following lines at the end of your configuration file: 

    from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
    process.muonDebugger =cms.EDAnalyzer("MuonHLTDebugger",
                                     MuonServiceProxy,
                                     triggerResults   = cms.untracked.InputTag("TriggerResults::TEST"),
                                     triggerSummary   = cms.untracked.InputTag("hltTriggerSummaryAOD::TEST"),
                                     L3Candidates     = cms.untracked.InputTag("hltIterL3OIL3MuonCandidates"),
                                     L2Candidates     = cms.untracked.InputTag("hltL2MuonCandidates"),
                                     L1Candidates     = cms.untracked.InputTag("hltGmtStage2Digis", "Muon"),
                                     MuonLinksTag     = cms.untracked.InputTag("hltIterL3MuonsFromL2LinksCombination","","TEST"),
                                     genParticlesTag  = cms.untracked.InputTag("genParticles"),
                                     muonTag          = cms.untracked.InputTag("muons"),
                                     triggerProcess   = cms.string("TEST"),
                                     triggerName      = cms.string("HLT_Mu50_v5"),
                                     theTracksOISeeds = cms.untracked.InputTag("hltIterL3OISeedsFromL2Muons","","TEST"),
                                     theTracksOICand  = cms.untracked.InputTag("hltIterL3OITrackCandidates","","TEST"),
                                     theTracksOINoHP  = cms.untracked.InputTag("hltIterL3OIMuCtfWithMaterialTracks","","TEST"),
                                     theTracksOI      = cms.untracked.InputTag("hltIterL3OIMuonTrackSelectionHighPurity","","TEST"),
                                     l1filterLabel    = cms.string("hltL1fL1sMu22Or25L1Filtered0"),
                                     l2filterLabel    = cms.string("hltL2fL1sMu22or25L1f0L2Filtered10Q"),
                                     l3filterLabel    = cms.string("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q"),
                                     debuglevel       = cms.untracked.uint32(2),
                                     runSharedHits    = cms.untracked.bool(False),
                                     isMC             = cms.untracked.bool(True)
    )
    process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonDebugger_MC_IterL3.root"),
                                   closeFileFast = cms.untracked.bool(False)
    )
    process.HLTValidation = cms.EndPath(
        process.muonDebugger
    )

## Plotter and other Tools
