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

First thing to do is to download a configuration file with the desired HLT paths, extract an HLT configuration from ConfDB using the hltGetConfiguration script, here an example: 

hltGetConfiguration /online/collisions/2017/2e34/v4.2/HLT/V6 --input file:/eos/cms/store/user/folguera/ROOT/BA7BEDED-948E-E711-AFC7-02163E011FAF.root  --paths HLTriggerFirstPath,HLT_IsoMu27_v14,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v13,HLT_Mu50_v12,HLT_DoubleMu43NoFiltersNoVtx_v3,HLT_Dimuon25_Jpsi_v13,DST_DoubleMu3_noVtx_CaloScouting_v5,HLTriggerFinalPath --output none --full --offline --data --unprescale --process TEST --globaltag auto:run2_hlt_GRun > hltConfig.py

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

Then you can run the configuration file with cmsRun locally, using the batch (Tools/submitJobs.py) or using crab (i.e. python Tools/crab_HLTMC_eff_multicrab_cfg.py). 

## Run debugger (for event-by-event comparison) 

You need to add the following lines at the end of your configuration file: 

from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
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

process.HLTValidation = cms.EndPath(
    process.muonDebugger
)

## Plotter and other Tols

The 
