import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.source = cms.Source("PoolSource",
                    fileNames = cms.untracked.vstring(
                      '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/FA0B284A-856B-E511-9378-02163E0133E8.root',
                      '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/FAC77F63-8F6B-E511-9911-02163E0145FD.root',
                      '/store/data/Run2015D/SingleMuon/RAW-RECO/ZMu-PromptReco-v3/000/258/158/00000/FE6D8959-7C6B-E511-B41B-02163E0137BC.root',
                    ),
                    secondaryFileNames = cms.untracked.vstring(),
#                     lumisToProcess = cms.untracked.VLuminosityBlockRange('258158:1-258158:1786'),

)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


process.muonNtuples =cms.EDAnalyzer("MuonNtuples",
                       offlineVtx               = cms.InputTag("offlinePrimaryVertices"),
                       offlineMuons             = cms.InputTag("muons"),
                       
                       triggerResult            = cms.untracked.InputTag("TriggerResults::HLT"),
                       triggerSummary           = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                       tagTriggerResult         = cms.untracked.InputTag("TriggerResults::HLT"),
                       tagTriggerSummary        = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                       
                       L3Candidates             = cms.untracked.InputTag("hltL3MuonCandidates"), 
                       L2Candidates             = cms.untracked.InputTag("hltL2MuonCandidates"), 
                       L1Candidates             = cms.untracked.InputTag("hltGmtStage2Digis", "Muon"), 
                       TkMuCandidates           = cms.untracked.InputTag("hltHighPtTkMuonCands"), 
                       NeutralDeposit           = cms.untracked.InputTag("hltMuonHcalPFClusterIsoForMuons"), 
                       PhotonsDeposit           = cms.untracked.InputTag("hltMuonEcalPFClusterIsoForMuons"), 
                       NeutralDeposit05         = cms.untracked.InputTag("hltMuonHcalPFClusterIsoForMuonsNoEffAreaVeto0p05"), 
                       PhotonsDeposit05         = cms.untracked.InputTag("hltMuonEcalPFClusterIsoForMuonsNoEffAreaVeto0p05"), 
                       NeutralDeposit1          = cms.untracked.InputTag("hltMuonHcalPFClusterIsoForMuonsNoEffAreaVeto0p1"), 
                       PhotonsDeposit1          = cms.untracked.InputTag("hltMuonEcalPFClusterIsoForMuonsNoEffAreaVeto0p1"), 
                       ChargedDeposit           = cms.untracked.InputTag("hltMuonTkRelIsolationCut0p09Map", "trkIsoDeposits", "HLT"), 
                       
                       RhoCorrectionOnline      = cms.untracked.InputTag("hltFixedGridRhoFastjetAllCaloForMuons"), # for now, same for tag and probe muons
                       RhoCorrectionOffline     = cms.untracked.InputTag("fixedGridRhoFastjetAllCalo"), 
                                              
                       offlineECalPFIso03       = cms.untracked.InputTag("muonEcalPFClusterIsolationProducer03"), 
                       offlineHCalPFIso03       = cms.untracked.InputTag("muonHcalPFClusterIsolationProducer03"), 
                       offlineECalPFIso04       = cms.untracked.InputTag("muonEcalPFClusterIsolationProducer04"), 
                       offlineHCalPFIso04       = cms.untracked.InputTag("muonHcalPFClusterIsolationProducer04"), 
                       
                       lumiScalerTag            = cms.untracked.InputTag("scalersRawToDigi"),
                       puInfoTag                = cms.untracked.InputTag("addPileupInfo"),

                       genParticlesTag          = cms.untracked.InputTag("genParticles"),
                       doOffline                = cms.untracked.bool(True)
                       )   
                       
process.mypath  = cms.Path(process.muonNtuples)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))   


process.MessageLogger = cms.Service("MessageLogger",
   destinations   = cms.untracked.vstring('cerr'),
   cerr           = cms.untracked.PSet(
       threshold      = cms.untracked.string('ERROR'),
   ),
#    debugModules  = cms.untracked.vstring('*')
)

