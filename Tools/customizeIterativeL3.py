import FWCore.ParameterSet.Config as cms


def customizeForOIWithSeedCleaningAsCascade(process):
    process.hltIterL3OITrackCandidates.TrajectoryCleaner = cms.string("hltESPTrajectoryCleanerBySharedHits")
    
    return process

def customizeForOIForGroupedCkfTrajectory(process):
    process.HLTPSetCkfTrajectoryFilterIterL3OI = cms.PSet( 
	  minimumNumberOfHits = cms.int32( 5 ),
	  ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
	  seedExtension = cms.int32( 0 ),
	  chargeSignificance = cms.double( -1.0 ),
	  pixelSeedExtension = cms.bool( False ),
	  strictSeedExtension = cms.bool( False ),
	  nSigmaMinPt = cms.double( 5.0 ),
	  maxCCCLostHits = cms.int32( 9999 ),
	  minPt = cms.double( 3.0 ),
	  maxConsecLostHits = cms.int32( 1 ),
	  extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
	  constantValueForLostHitsFractionFilter = cms.double( 10.0 ),
	  seedPairPenalty = cms.int32( 0 ),
	  maxNumberOfHits = cms.int32( -1 ),
	  minNumberOfHitsForLoopers = cms.int32( 13 ),
	  minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
	  minNumberOfHitsPerLoop = cms.int32( 4 ),
	  minHitsMinPt = cms.int32( 3 ),
	  maxLostHitsFraction = cms.double( 999.0 ),
	  maxLostHits = cms.int32( 1 )
          )
    process.HLTPSetGroupedCkfTrajectoryBuilderIterL3ForOI = cms.PSet( 
	  rescaleErrorIfFail = cms.double( 1.0 ),
	  keepOriginalIfRebuildFails = cms.bool( False ),
	  lockHits = cms.bool( True ),
	  propagatorOpposite = cms.string( "PropagatorWithMaterialOpposite" ),
	  trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfTrajectoryFilterIterL3OI" ) ),
	  maxCand = cms.int32( 5 ),
	  estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator30" ),
	  intermediateCleaning = cms.bool( True ),
	  bestHitOnly = cms.bool( True ),
	  deltaEta = cms.double( -1.0 ),
	  useSeedLayer = cms.bool( False ),
	  useSameTrajFilter = cms.bool( True ),
	  MeasurementTrackerName = cms.string( "hltSiStripClusters" ),
	  ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
	  lostHitPenalty = cms.double( 30.0 ),
	  requireSeedHitsInRebuild = cms.bool( False ),
	  TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
	  propagatorAlong = cms.string( "PropagatorWithMaterial" ),
	  minNrOfHitsForRebuild = cms.int32( 5 ),
	  alwaysUseInvalidHits = cms.bool( True ),
	  inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTPSetCkfTrajectoryFilterIterL3OI" ) ),
	  foundHitBonus = cms.double( 1000.0 ),
	  propagatorProximity = cms.string( "SteppingHelixPropagatorAny" ),
	  updator = cms.string( "hltESPKFUpdator" ),
	  deltaPhi = cms.double( -1.0 )
          )
    process.hltIterL3OITrackCandidates.TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTPSetGroupedCkfTrajectoryBuilderIterL3ForOI" ) )
    process.hltIterL3OITrackCandidates.TrajectoryBuilder = cms.string( "GroupedCkfTrajectoryBuilder" )

    return process

def customizeForOIWithHitBasedOnly(process):
    process.hltIterL3OISeedsFromL2Muons.UseHitLessSeeds = cms.bool( False )
    return process

def customizeTSGForOIHitsToTry1(process):
    process.hltIterL3OISeedsFromL2Muons.hitsToTry = cms.int32( 1 )
    process.hltIterL3OISeedsFromL2Muons.layersToTry = cms.int32( 5 )
    process.hltIterL3OISeedsFromL2Muons.maxSeeds = cms.uint32( 5 )
    return process

def customizeTSGForOIHitsToTry2Layer5(process):
    process.hltIterL3OISeedsFromL2Muons.hitsToTry = cms.int32( 2 )
    process.hltIterL3OISeedsFromL2Muons.layersToTry = cms.int32( 5 )
    process.hltIterL3OISeedsFromL2Muons.maxSeeds = cms.uint32( 10 )
    return process

def customizeTSGForOISmallerOverlap(process):
    process.hltIterL3OISeedsFromL2Muons.maxEtaForTOB = cms.double( 1.4 )
    process.hltIterL3OISeedsFromL2Muons.minEtaForTEC = cms.double( 0.9 )
    return process
    
##def customizeForOIWithoutSeedCleaner(process):
##    process.hltIterL3OITrackCandidates.RedundantSeedCleaner = cms.string( "none" ) 
##    return process

def customizeForOIWithOfflineParametersForTrackCandidate(process):
    process.hltIterL3OITrackCandidates.numHitsForSeedCleaner = cms.int32(50)
    process.hltIterL3OITrackCandidates.onlyPixelHitsForSeedCleaner = cms.bool(False)
    
    return process


def customizeTSGForOIWithPropagator(process):
	process.hltIterL3OISeedsFromL2Muons = cms.EDProducer( "TSGForOI",
	    hitsToTry = cms.int32( 3 ),
	    adjustErrorsDynamicallyForHitless = cms.bool( True ),
	    SF4 = cms.double( 7.0 ),
	    SF5 = cms.double( 10.0 ),
	    SF2 = cms.double( 4.0 ),
	    SF3 = cms.double( 5.0 ),
	    SF1 = cms.double( 3.0 ),
	    minEtaForTEC = cms.double( 0.7 ),
	    fixedErrorRescaleFactorForHits = cms.double( 3.0 ),
	    maxSeeds = cms.uint32( 5 ),
	    maxEtaForTOB = cms.double( 1.8 ),
	    pT3 = cms.double( 70.0 ),
	    pT2 = cms.double( 30.0 ),
	    pT1 = cms.double( 13.0 ),
	    layersToTry = cms.int32( 2 ),
	    fixedErrorRescaleFactorForHitless = cms.double( 10.0 ),
	    MeasurementTrackerEvent = cms.InputTag( "hltSiStripClusters" ),
	    adjustErrorsDynamicallyForHits = cms.bool( True ),
	    src = cms.InputTag( 'hltL2Muons','UpdatedAtVtx' ),
	    tsosDiff = cms.double( 0.03 ),
	    eta1 = cms.double( 1.0 ),
	    eta2 = cms.double( 1.4 ),
	    UseHitLessSeeds = cms.bool( True ),
	    UseStereoLayersInTEC = cms.bool( False ),
	    estimator = cms.string( "hltESPChi2MeasurementEstimator100" ),
	    debug = cms.untracked.bool( False ),
            propagatorAlongName = cms.string( "PropagatorWithMaterial" ),
            propagatorOppositeName = cms.string( "PropagatorWithMaterial" ),
            propagatorHitlessName = cms.string( "PropagatorWithMaterial" ),
            propagatorHitBasedName = cms.string( "PropagatorWithMaterial" ),
	)
    

def customizeForOILooseFilter(process):
	process.hltIterL3OIMuonTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
	    src = cms.InputTag( "hltIterL3OIMuCtfWithMaterialTracks" ),
	    GBRForestLabel = cms.string( "" ),
	    beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
	    vertices = cms.InputTag( "Notused" ),
	    qualityCuts = cms.vdouble( -0.7, 0.1, 0.4 ),
	    mva = cms.PSet( 
	      minPixelHits = cms.vint32( 0, 0, 1 ),
	      maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 100.0 ),
	      dr_par = cms.PSet( 
	        d0err = cms.vdouble( 0.003, 0.003, 3.40282346639E38 ),
	        dr_par2 = cms.vdouble( 0.3, 0.3, 3.40282346639E38 ),
	        dr_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
	        dr_exp = cms.vint32( 4, 4, 2147483647 ),
	        d0err_par = cms.vdouble( 0.001, 0.001, 3.40282346639E38 )
	      ),
	      maxLostLayers = cms.vint32( 4, 3, 4 ),
	      min3DLayers = cms.vint32( 1, 2, 1 ),
	      dz_par = cms.PSet( 
	        dz_par1 = cms.vdouble( 0.4, 0.4, 3.40282346639E38 ),
	        dz_par2 = cms.vdouble( 0.35, 0.35, 3.40282346639E38 ),
	        dz_exp = cms.vint32( 4, 4, 2147483647 )
	      ),
	      minNVtxTrk = cms.int32( 2 ),
	      maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
	      minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
	      maxChi2 = cms.vdouble( 3.40282346639E38, 3.40282346639E38, 3.40282346639E38 ),
	      maxChi2n = cms.vdouble( 10.0, 1.0, 1.0 ),
	      maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
	      minLayers = cms.vint32( 3, 5, 3 )
	    ),
	    ignoreVertices = cms.bool( True ),
	    GBRForestFileName = cms.string( "" )
	)


    


