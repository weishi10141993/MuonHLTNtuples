import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:myfile.root'
    )
)

process.muonDebugger =cms.EDAnalyzer("MuonHLTDebugger",
                                    triggerResults  = cms.untracked.InputTag("TriggerResults::SFHLT"),
                                    triggerSummary  = cms.untracked.InputTag("hltTriggerSummaryAOD::SFHLT"),
                                    L3Candidates    = cms.untracked.InputTag("hltL3MuonCandidates"), 
                                    L2Candidates    = cms.untracked.InputTag("hltL2MuonCandidates"), 
                                    L1Candidates    = cms.untracked.InputTag("hltGmtStage2Digis", "Muon"), 
                                    genParticlesTag = cms.untracked.InputTag("genParticles"),
                                    triggerProcess  = cms.untracker.InputTag("SFHLT"),
                                    triggerName     = cms.untracker.InputTag("HLT_Mu50_v4"),
                                    l1filterLabel   = cms.untracker.InputTag("hltL1fL1sMu22Or25L1Filtered0"),
                                    l2filterLabel   = cms.untracker.InputTag("hltL2fL1sMu22Or25L1f0L2Filtered10Q"),
                                    l3filterLabel   = cms.untracker.InputTag("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q")
                                    )   

process.mypath  = cms.Path(process.muonDebugger)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("muonNtuple.root"),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))   


