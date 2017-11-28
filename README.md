# MuonHLTNtuples

cmsrel CMSSW_X_Y_Z  (-> update to the appropriate CMSSW version!)   
cd CMSSW_X_Y_Z/src  
cmsenv    
git cms-addpkg HLTrigger/Configuration    
git clone git@github.com:folguera/MuonHLTNtuples.git    
cd MuonHLTNtuples/  
git checkout -b YOUR_BRANCH_NAME    
cd ..  
scramv1 b   

First thing to do is to download a configuration file with the desired HLT paths, here an example: 

hltGetConfiguration /online/collisions/2017/2e34/v3.0/HLT/V16 --input file:/path/to/file/with/RAW/info.root  --paths HLTriggerFirstPath,HLT_IsoMu27_v13,HLT_Mu50_v*,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v12,HLT_DoubleMu43NoFiltersNoVtx_v*,HLT_Dimuon25_Jpsi_v*,DST_DoubleMu3_noVtx_CaloScouting_v*,HLTriggerFinalPath --output none --full --offline --data --unprescale --process TEST --globaltag auto:run2_hlt_GRun

Then at the end of it you should add different lines to run the Ntuple-producer or the Debugger: 

## Produce Ntuples: 
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

Then you can run the configuration file with cmsRun locally, using the bacth (Tools/submitJobs.py) or using crab (i.e. Tools/crab_HLTMC_eff_multicrab_cfg.py). 

## Run debugger (for event-by-event comparison) 

You need to add the following lines at the end of your configuration file: 

##from RecoMuon.TrackingTools.MuonServiceProxy_cff import *

process.muonDebugger =cms.EDAnalyzer("MuonHLTDebugger",
                                     MuonServiceProxy,
                                     triggerResults  = cms.untracked.InputTag("TriggerResults::SFHLT"),
                                     triggerSummary  = cms.untracked.InputTag("hltTriggerSummaryAOD::SFHLT"),
                                     L3Candidates    = cms.untracked.InputTag("hltNewL3MuonCandidates"),
                                     L2Candidates    = cms.untracked.InputTag("hltL2MuonCandidates"),
                                     L1Candidates    = cms.untracked.InputTag("hltGmtStage2Digis", "Muon"),
                                     MuonLinksTag    = cms.untracked.InputTag("hltNewL3MuonsLinksCombination","","SFHLT"),
                                     genParticlesTag = cms.untracked.InputTag("genParticles"),
                                     muonTag         = cms.untracked.InputTag("muons"),
                                     triggerProcess  = cms.string("SFHLT"),
                                     triggerName     = cms.string("HLT_Mu50_v5"),
                                     l1filterLabel   = cms.string("hltL1fL1sMu22Or25L1Filtered0"),
                                     l2filterLabel   = cms.string("hltL2fL1sMu22Or25L1f0L2Filtered10Q"),
                                     l3filterLabel   = cms.string("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"),
                                     debuglevel      = cms.untracked.uint32(0),
                                     isMC            = cms.untracked.bool(True)
                                     )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonDebugger_MC_IterL3.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

## Plotter and other Tols

The 
