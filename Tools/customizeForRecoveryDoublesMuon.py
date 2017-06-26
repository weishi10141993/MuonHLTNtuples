
import FWCore.ParameterSet.Config as cms

def customizeForPixelTriplets(process):
    process.hltPixelTripletsClustersRefRemoval = cms.EDProducer( "TrackClusterRemover",
    	trackClassifier = cms.InputTag( '','QualityMasks' ),
    	minNumberOfLayersWithMeasBeforeFiltering = cms.int32( 0 ),
    	maxChi2 = cms.double( 3000.0 ),
    	trajectories = cms.InputTag( "hltPixelTracks" ),
    	oldClusterRemovalInfo = cms.InputTag( "" ),
    	stripClusters = cms.InputTag( "hltSiStripRawToClustersFacility" ),
    	#stripClusters = cms.InputTag( "none" ),
    	overrideTrkQuals = cms.InputTag( "" ),
    	pixelClusters = cms.InputTag( "hltSiPixelClusters" ),
    	#TrackQuality = cms.string( "highPurity" )
    	TrackQuality = cms.string( "undefQuality" )
    )
    process.hltPixelLayerTriplets = cms.EDProducer( "SeedingLayersEDProducer",
    	layerList = cms.vstring(
       # 	'BPix1+BPix2+BPix3',
       # 	'BPix2+BPix3+BPix4',
       # 	'BPix1+BPix3+BPix4',
       # 	'BPix1+BPix2+BPix4',
       # 	'BPix2+BPix3+FPix1_pos',
       # 	'BPix2+BPix3+FPix1_neg',
       # 	'BPix1+BPix2+FPix1_pos',
       # 	'BPix1+BPix2+FPix1_neg',
       # 	'BPix2+FPix1_pos+FPix2_pos',
       # 	'BPix2+FPix1_neg+FPix2_neg',
       # 	'BPix1+FPix1_pos+FPix2_pos',
       # 	'BPix1+FPix1_neg+FPix2_neg',
       # 	'FPix1_pos+FPix2_pos+FPix3_pos',
       # 	'FPix1_neg+FPix2_neg+FPix3_neg'
	      	'BPix1+BPix2+BPix3',
		'BPix2+BPix3+BPix4',
      		'BPix1+BPix3+BPix4',
      		'BPix1+BPix2+BPix4',
      		'BPix2+BPix3+FPix1_pos',
      		'BPix2+BPix3+FPix1_neg',
      		'BPix1+BPix2+FPix1_pos',
      		'BPix1+BPix2+FPix1_neg',
      		'BPix2+FPix1_pos+FPix2_pos',
      		'BPix2+FPix1_neg+FPix2_neg',
      		'BPix1+FPix1_pos+FPix2_pos',
      		'BPix1+FPix1_neg+FPix2_neg',
      		'FPix1_pos+FPix2_pos+FPix3_pos',
      		'FPix1_neg+FPix2_neg+FPix3_neg',
		'BPix1+BPix3+FPix1_pos',
		'BPix1+BPix2+FPix2_pos',	
		'BPix1+BPix3+FPix1_neg',
		'BPix1+BPix2+FPix2_neg',
		'BPix1+FPix2_neg+FPix3_neg',
		'BPix1+FPix1_neg+FPix3_neg',  
		'BPix1+FPix2_pos+FPix3_pos',
		'BPix1+FPix1_pos+FPix3_pos',  

	),		   
    	MTOB = cms.PSet(  ),
    	TEC = cms.PSet(  ),
    	MTID = cms.PSet(  ),
    	FPix = cms.PSet( 
      		hitErrorRPhi = cms.double( 0.0051 ),
      		TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      		useErrorsFromParam = cms.bool( True ),
      		hitErrorRZ = cms.double( 0.0036 ),
      		HitProducer = cms.string( "hltSiPixelRecHits" ),
      		skipClusters = cms.InputTag( "hltPixelTripletsClustersRefRemoval" )
    	),
    	MTEC = cms.PSet(  ),
    	MTIB = cms.PSet(  ),
    	TID = cms.PSet(  ),
    	TOB = cms.PSet(  ),
    	BPix = cms.PSet( 
      		hitErrorRPhi = cms.double( 0.0027 ),
      		TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      		useErrorsFromParam = cms.bool( True ),
      		hitErrorRZ = cms.double( 0.006 ),
      		HitProducer = cms.string( "hltSiPixelRecHits" ),
      		skipClusters = cms.InputTag( "hltPixelTripletsClustersRefRemoval" )
    	),
    	TIB = cms.PSet(  )
    )
    
    process.hltPixelTracksTrackingRegionsForTriplets = cms.EDProducer( "PointSeededTrackingRegionsEDProducer",
    	RegionPSet = cms.PSet( 
     		vertexCollection = cms.InputTag( "none" ),
    		zErrorVetex = cms.double( 0.1 ),
      		beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      		zErrorBeamSpot = cms.double( 15.0 ),
      		maxNVertices = cms.int32( 10 ),
      		maxNRegions = cms.int32( 100 ),
      		nSigmaZVertex = cms.double( 3.0 ),
      		nSigmaZBeamSpot = cms.double( 3.0 ),
      		ptMin = cms.double( 0.8 ),
      		mode = cms.string( "BeamSpotFixed" ),
      		searchOpt = cms.bool( False ),
      		whereToUseMeasurementTracker = cms.string( "never" ),
      		originRadius = cms.double( 0.1 ),
      		measurementTrackerName = cms.InputTag( "hltIter4MaskedMeasurementTrackerEvent" ),
      		precise = cms.bool( True ),
      		deltaEta = cms.double( 1.2 ),
      		deltaPhi = cms.double( 0.5 ),
		points = cms.VPSet(
    			cms.PSet(
				eta = cms.double(0.4),
				phi = cms.double(3.0)
	
			),
			cms.PSet(
				eta = cms.double(1.8),
				phi = cms.double(2.0)
	
			),
  			cms.PSet(
				eta = cms.double(-1.8),
				phi = cms.double(0.25)
	
			),

		)
    	)

    )



    process.hltPixelTracksHitDoubletsForTriplets = cms.EDProducer( "HitPairEDProducer",
    	trackingRegions = cms.InputTag( "hltPixelTracksTrackingRegionsForTriplets" ),
    	layerPairs = cms.vuint32( 0 ),
    	clusterCheck = cms.InputTag( "" ),
    	produceSeedingHitSets = cms.bool( False ),
    	produceIntermediateHitDoublets = cms.bool( True ),
    	maxElement = cms.uint32( 0 ),
    	seedingLayers = cms.InputTag( "hltPixelLayerTriplets" )
    )
    process.hltPixelTracksHitTriplets = cms.EDProducer( "CAHitTripletEDProducer",
    	CAThetaCut = cms.double( 0.002 ),
    	SeedComparitorPSet = cms.PSet( 
      	clusterShapeHitFilter = cms.string( "ClusterShapeHitFilter" ),
      	ComponentName = cms.string( "LowPtClusterShapeSeedComparitor" ),
      	clusterShapeCacheSrc = cms.InputTag( "hltSiPixelClustersCache" )
    	),
    	extraHitRPhitolerance = cms.double( 0.032 ),
    	doublets = cms.InputTag( "hltPixelTracksHitDoubletsForTriplets" ),
    	CAHardPtCut = cms.double( 0.0 ),
    		maxChi2 = cms.PSet( 
      		value2 = cms.double( 50.0 ),
      		value1 = cms.double( 200.0 ),
      		pt1 = cms.double( 0.7 ),
      		enabled = cms.bool( False ),
      		pt2 = cms.double( 2 )
    	),
    	CAPhiCut = cms.double( 0.2 ),
    	useBendingCorrection = cms.bool( True ),
    )


    process.hltPixelTracksFromTriplets = cms.EDProducer( "PixelTrackProducer",
    	Filter = cms.InputTag( "hltPixelTracksFilter" ),
    	Cleaner = cms.string( "hltPixelTracksCleanerBySharedHits" ),
    	passLabel = cms.string( "" ),
    	Fitter = cms.InputTag( "hltPixelTracksFitter" ),
    	SeedingHitSets = cms.InputTag( "hltPixelTracksHitTriplets" )
    )

    process.hltPixelTrackMerged = cms.EDProducer( "TrackListMerger",
    	ShareFrac = cms.double( 0.19 ),
    	writeOnlyTrkQuals = cms.bool( False ),
    	MinPT = cms.double( 0.05 ),
    	allowFirstHitShare = cms.bool( True ),
    	copyExtras = cms.untracked.bool( True ),
    	Epsilon = cms.double( -0.001 ),
    	selectedTrackQuals = cms.VInputTag( 'hltPixelTracks','hltPixelTracksFromTriplets' ),
    	indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    	MaxNormalizedChisq = cms.double( 1000.0 ),
    	copyMVA = cms.bool( False ),
    	FoundHitBonus = cms.double( 5.0 ),
    	setsToMerge = cms.VPSet( 
      		cms.PSet(  pQual = cms.bool( False ),
        		tLists = cms.vint32( 0, 1 )
      		)
    	),
    	MinFound = cms.int32( 3 ),
    	hasSelector = cms.vint32( 0, 0 ),
    	TrackProducers = cms.VInputTag( 'hltPixelTracks','hltPixelTracksFromTriplets' ),
    	LostHitPenalty = cms.double( 20.0 ),
    	trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    	newQuality = cms.string( "confirmed" )
    )

   # process.hltIter0PFLowPixelSeedsFromPixelTracks.InputCollection = cms.InputTag('hltPixelTrackMerged')
   
    process.hltPixelVertices.TrackCollection = cms.InputTag( "hltPixelTrackMerged" )
 
    process.HLTRecopixelvertexingSequence = cms.Sequence(process.hltPixelTracksFilter+process.hltPixelTracksFitter+process.hltPixelTracksTrackingRegions+process.hltPixelLayerQuadruplets+process.hltPixelTracksHitDoublets+process.hltPixelTracksHitQuadruplets+process.hltPixelTracks+process.hltPixelTripletsClustersRefRemoval+process.hltPixelTracksTrackingRegionsForTriplets+process.hltPixelLayerTriplets+process.hltPixelTracksHitDoubletsForTriplets+process.hltPixelTracksHitTriplets+process.hltPixelTracksFromTriplets+process.hltPixelTrackMerged+process.hltPixelVertices+process.hltTrimmedPixelVertices)
    process.MC_ReducedIterativeTracking_v3 = cms.Path( process.HLTBeginSequence + process.hltPreMCReducedIterativeTracking + process.HLTRecoJetSequenceAK4PrePF + process.HLTDoLocalPixelSequence + process.HLTDoLocalStripSequence + process.HLTRecopixelvertexingSequence + process.HLTIterativeTrackingIter02 + process.HLTEndSequence )
    #pr
    return process

def customizeForRecoveryIter2(process):
    ### Add layer combinations per recommendation from Matti
    process.hltIter2PixelLayerTriplets.layerList = cms.vstring( 
	      	'BPix1+BPix2+BPix3',
		'BPix2+BPix3+BPix4',
      		'BPix1+BPix3+BPix4',
      		'BPix1+BPix2+BPix4',
      		'BPix2+BPix3+FPix1_pos',
      		'BPix2+BPix3+FPix1_neg',
      		'BPix1+BPix2+FPix1_pos',
      		'BPix1+BPix2+FPix1_neg',
      		'BPix2+FPix1_pos+FPix2_pos',
      		'BPix2+FPix1_neg+FPix2_neg',
      		'BPix1+FPix1_pos+FPix2_pos',
      		'BPix1+FPix1_neg+FPix2_neg',
      		'FPix1_pos+FPix2_pos+FPix3_pos',
      		'FPix1_neg+FPix2_neg+FPix3_neg',
		'BPix1+BPix3+FPix1_pos',
		'BPix1+BPix2+FPix2_pos',	
		'BPix1+BPix3+FPix1_neg',
		'BPix1+BPix2+FPix2_neg',
		'BPix1+FPix2_neg+FPix3_neg',
		'BPix1+FPix1_neg+FPix3_neg',  
		'BPix1+FPix2_pos+FPix3_pos',
		'BPix1+FPix1_pos+FPix3_pos',  
    )



    return process
def customizeForRecoveryTriplets(process):

    process.hltIter3ClustersRefRemoval = cms.EDProducer( "TrackClusterRemover",
    	trackClassifier = cms.InputTag( '','QualityMasks' ),
    	minNumberOfLayersWithMeasBeforeFiltering = cms.int32( 0 ),
    	maxChi2 = cms.double( 16.0 ),
    	trajectories = cms.InputTag( "hltIter2PFlowTrackSelectionHighPurity" ),
    	oldClusterRemovalInfo = cms.InputTag( "hltIter2ClustersRefRemoval" ),
    	stripClusters = cms.InputTag( "hltSiStripRawToClustersFacility" ),
    	overrideTrkQuals = cms.InputTag( "" ),
    	pixelClusters = cms.InputTag( "hltSiPixelClusters" ),
    	TrackQuality = cms.string( "highPurity" )
    )
    process.hltIter3MaskedMeasurementTrackerEvent = cms.EDProducer( "MaskedMeasurementTrackerEventProducer",
    	clustersToSkip = cms.InputTag( "hltIter3ClustersRefRemoval" ),
    	OnDemand = cms.bool( False ),
    	src = cms.InputTag( "hltSiStripClusters" )
    )

    process.hltIter3PixelLayerTriplets = cms.EDProducer( "SeedingLayersEDProducer",
    	layerList = cms.vstring( 'BPix1+BPix2+BPix3',
      'BPix2+BPix3+BPix4',
      'BPix1+BPix3+BPix4',
      'BPix1+BPix2+BPix4',
      'BPix2+BPix3+FPix1_pos',
      'BPix2+BPix3+FPix1_neg',
      'BPix1+BPix2+FPix1_pos',
      'BPix1+BPix2+FPix1_neg',
      'BPix2+FPix1_pos+FPix2_pos',
      'BPix2+FPix1_neg+FPix2_neg',
      'BPix1+FPix1_pos+FPix2_pos',
      'BPix1+FPix1_neg+FPix2_neg',
      'FPix1_pos+FPix2_pos+FPix3_pos',
      'FPix1_neg+FPix2_neg+FPix3_neg' ),
    	MTOB = cms.PSet(  ),
    	TEC = cms.PSet(  ),
    	MTID = cms.PSet(  ),
    	FPix = cms.PSet( 
      		hitErrorRPhi = cms.double( 0.0051 ),
      		TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      		skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
      		useErrorsFromParam = cms.bool( True ),
      		hitErrorRZ = cms.double( 0.0036 ),
      		HitProducer = cms.string( "hltSiPixelRecHits" )
    	),
    	MTEC = cms.PSet(  ),
    	MTIB = cms.PSet(  ),
    	TID = cms.PSet(  ),
    	TOB = cms.PSet(  ),
    	BPix = cms.PSet( 
      		hitErrorRPhi = cms.double( 0.0027 ),
      		TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
      		skipClusters = cms.InputTag( "hltIter2ClustersRefRemoval" ),
      		useErrorsFromParam = cms.bool( True ),
      		hitErrorRZ = cms.double( 0.006 ),
      		HitProducer = cms.string( "hltSiPixelRecHits" )
    	),
    	TIB = cms.PSet(  )
    )

    process.hltIter3PixelLayerTriplets.layerList = cms.vstring( 
	      	'BPix1+BPix2+BPix3',
		'BPix2+BPix3+BPix4',
      		'BPix1+BPix3+BPix4',
      		'BPix1+BPix2+BPix4',
      		'BPix2+BPix3+FPix1_pos',
      		'BPix2+BPix3+FPix1_neg',
      		'BPix1+BPix2+FPix1_pos',
      		'BPix1+BPix2+FPix1_neg',
      		'BPix2+FPix1_pos+FPix2_pos',
      		'BPix2+FPix1_neg+FPix2_neg',
      		'BPix1+FPix1_pos+FPix2_pos',
      		'BPix1+FPix1_neg+FPix2_neg',
      		'FPix1_pos+FPix2_pos+FPix3_pos',
      		'FPix1_neg+FPix2_neg+FPix3_neg',
[		'BPix1+BPix3+FPix1_pos',
		'BPix1+BPix2+FPix2_pos',	
		'BPix1+BPix3+FPix1_neg',
		'BPix1+BPix2+FPix2_neg',
		'BPix1+FPix2_neg+FPix3_neg',
		'BPix1+FPix1_neg+FPix3_neg',  
		'BPix1+FPix2_pos+FPix3_pos',
		'BPix1+FPix1_pos+FPix3_pos',  
    )

    #process.hltIter2PFlowPixelTrackingRegions = cms.EDProducer( "GlobalTrackingRegionFromBeamSpotEDProducer",
    #	RegionPSet = cms.PSet(
     # 		nSigmaZ = cms.double( 0.0 ),
      #		beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      #		ptMin = cms.double( 0.9 ),
      #		originHalfLength = cms.double( 24.0 ),
      #		originRadius = cms.double( 0.2 ),
      #		precise = cms.bool( True ),
      #		useMultipleScattering = cms.bool( False )
   #	)
   # )
    process.hltIter3PFlowPixelTrackingRegions = cms.EDProducer( "PointSeededTrackingRegionsEDProducer",
    	RegionPSet = cms.PSet( 
     		vertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
    		zErrorVetex = cms.double( 0.1 ),
      		beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      		zErrorBeamSpot = cms.double( 15.0 ),
      		maxNVertices = cms.int32( 10 ),
      		maxNRegions = cms.int32( 100 ),
      		nSigmaZVertex = cms.double( 3.0 ),
      		nSigmaZBeamSpot = cms.double( 3.0 ),
      		ptMin = cms.double( 0.8 ),
      		mode = cms.string( "VerticesFixed" ),
      		searchOpt = cms.bool( False ),
      		whereToUseMeasurementTracker = cms.string( "ForSiStrips" ),
      		originRadius = cms.double( 0.05 ),
      		measurementTrackerName = cms.InputTag( "hltIter3MaskedMeasurementTrackerEvent" ),
      		precise = cms.bool( True ),
      		deltaEta = cms.double( 0.8 ),
      		deltaPhi = cms.double( 0.6 ),
		points = cms.VPSet(
     			cms.PSet(
				eta = cms.double(1.8),
				phi = cms.double(2.0)
	
			),
  			cms.PSet(
				eta = cms.double(-1.8),
				phi = cms.double(0.25)
	
			),

		)
    	)

    )

    process.hltIter3PFlowPixelClusterCheck = cms.EDProducer( "ClusterCheckerEDProducer",
    	cut = cms.string( "" ),
    	silentClusterCheck = cms.untracked.bool( False ),
    	MaxNumberOfCosmicClusters = cms.uint32( 50000 ),
    	PixelClusterCollectionLabel = cms.InputTag( "hltSiPixelClusters" ),
    	doClusterCheck = cms.bool( False ),
    	MaxNumberOfPixelClusters = cms.uint32( 40000 ),
    	ClusterCollectionLabel = cms.InputTag( "hltSiStripClusters" )
    )
    process.hltIter3PFlowPixelHitDoublets = cms.EDProducer( "HitPairEDProducer",
    	trackingRegions = cms.InputTag( "hltIter3PFlowPixelTrackingRegions" ),
    	layerPairs = cms.vuint32( 0, 1),
    	clusterCheck = cms.InputTag( "hltIter3PFlowPixelClusterCheck" ),
    	produceSeedingHitSets = cms.bool( False ),
    	produceIntermediateHitDoublets = cms.bool( True ),
    	maxElement = cms.uint32( 0 ),
    	seedingLayers = cms.InputTag( "hltIter3PixelLayerTriplets" )
    )
    process.hltIter3PFlowPixelHitTriplets = cms.EDProducer( "CAHitTripletEDProducer",
    	CAHardPtCut = cms.double( 0.3 ),
    	SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
    	extraHitRPhitolerance = cms.double( 0.032 ),
    	doublets = cms.InputTag( "hltIter3PFlowPixelHitDoublets" ),
    	CAThetaCut = cms.double( 0.004 ),
    	maxChi2 = cms.PSet( 
      		value2 = cms.double( 50.0 ),
      		value1 = cms.double( 100.0 ),
      		pt1 = cms.double( 0.8 ),
      		enabled = cms.bool( True ),
      		pt2 = cms.double( 8.0 )
        ),
    	CAPhiCut = cms.double( 0.1 ),
    	useBendingCorrection = cms.bool( True )
    )
    process.hltIter3PFlowPixelSeeds = cms.EDProducer( "SeedCreatorFromRegionConsecutiveHitsTripletOnlyEDProducer",
    	SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
    	forceKinematicWithRegionDirection = cms.bool( False ),
    	magneticField = cms.string( "ParabolicMf" ),
    	SeedMomentumForBOFF = cms.double( 5.0 ),
    	OriginTransverseErrorMultiplier = cms.double( 1.0 ),
    	TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    	MinOneOverPtError = cms.double( 1.0 ),
    	seedingHitSets = cms.InputTag( "hltIter3PFlowPixelHitTriplets" ),
    	propagator = cms.string( "PropagatorWithMaterialParabolicMf" )
    )
    process.hltIter3PFlowCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
    	src = cms.InputTag( "hltIter3PFlowPixelSeeds" ),
    	maxSeedsBeforeCleaning = cms.uint32( 1000 ),
    	SimpleMagneticField = cms.string( "ParabolicMf" ),
    	TransientInitialStateEstimatorParameters = cms.PSet( 
      		propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      		numberMeasurementsForFit = cms.int32( 4 ),
      		propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    	),
    	TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    	MeasurementTrackerEvent = cms.InputTag( "hltIter3MaskedMeasurementTrackerEvent" ),
    	cleanTrajectoryAfterInOut = cms.bool( False ),
    	useHitsSplitting = cms.bool( False ),
    	RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
    	doSeedingRegionRebuilding = cms.bool( False ),
    	maxNSeeds = cms.uint32( 100000 ),
	produceSeedStopReasons = cms.bool( False ),
    	TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter2GroupedCkfTrajectoryBuilderIT" ) ),
    	NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    	TrajectoryBuilder = cms.string( "" )
    )
    process.hltIter3PFlowCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    	src = cms.InputTag( "hltIter3PFlowCkfTrackCandidates" ),
    	SimpleMagneticField = cms.string( "ParabolicMf" ),
    	clusterRemovalInfo = cms.InputTag( "" ),
    	beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    	MeasurementTrackerEvent = cms.InputTag( "hltIter3MaskedMeasurementTrackerEvent" ),
    	Fitter = cms.string( "hltESPFittingSmootherIT" ),
    	useHitsSplitting = cms.bool( False ),
    	MeasurementTracker = cms.string( "" ),
    	AlgorithmName = cms.string( "hltIter2" ),
    	alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    	NavigationSchool = cms.string( "" ),
    	TrajectoryInEvent = cms.bool( False ),
    	TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    	GeometricInnerState = cms.bool( True ),
    	useSimpleMF = cms.bool( True ),
    	Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" )
    )
    process.hltIter3PFlowTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    	src = cms.InputTag( "hltIter3PFlowCtfWithMaterialTracks" ),
    	GBRForestLabel = cms.string( "" ),
    	beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    	vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    	qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    	mva = cms.PSet( 
      	minPixelHits = cms.vint32( 0, 0, 0 ),
      	maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
      	dr_par = cms.PSet( 
        	d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
        	dr_par2 = cms.vdouble( 3.40282346639E38, 0.3, 0.3 ),
        	dr_par1 = cms.vdouble( 3.40282346639E38, 0.4, 0.4 ),
        	dr_exp = cms.vint32( 4, 4, 4 ),
        	d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
      	),
      	maxLostLayers = cms.vint32( 1, 1, 1 ),
      	min3DLayers = cms.vint32( 0, 0, 0 ),
      	dz_par = cms.PSet( 
        	dz_par1 = cms.vdouble( 3.40282346639E38, 0.4, 0.4 ),
        	dz_par2 = cms.vdouble( 3.40282346639E38, 0.35, 0.35 ),
        	dz_exp = cms.vint32( 4, 4, 4 )
      	),
      	minNVtxTrk = cms.int32( 3 ),
      	maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      	minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      	maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
      	maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      	maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      	minLayers = cms.vint32( 3, 3, 3 )
    	),
    	GBRForestFileName = cms.string( "" )
    )
    process.hltIter3PFlowTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
    	minQuality = cms.string( "highPurity" ),
        copyExtras = cms.untracked.bool( True ),
    	copyTrajectories = cms.untracked.bool( False ),
    	originalSource = cms.InputTag( "hltIter3PFlowCtfWithMaterialTracks" ),
    	originalQualVals = cms.InputTag( 'hltIter3PFlowTrackCutClassifier','QualityMasks' ),
    	originalMVAVals = cms.InputTag( 'hltIter3PFlowTrackCutClassifier','MVAValues' )
    )
    process.hltIter3Merged = cms.EDProducer( "TrackListMerger",
    	ShareFrac = cms.double( 0.19 ),
    	writeOnlyTrkQuals = cms.bool( False ),
    	MinPT = cms.double( 0.05 ),
    	allowFirstHitShare = cms.bool( True ),
    	copyExtras = cms.untracked.bool( True ),
    	Epsilon = cms.double( -0.001 ),
    	selectedTrackQuals = cms.VInputTag( 'hltIter2Merged','hltIter3PFlowTrackSelectionHighPurity' ),
    	indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    	MaxNormalizedChisq = cms.double( 1000.0 ),
    	copyMVA = cms.bool( False ),
    	FoundHitBonus = cms.double( 5.0 ),
   	 setsToMerge = cms.VPSet( 
      		cms.PSet(  pQual = cms.bool( False ),
        		tLists = cms.vint32( 0, 1 )
      		)
    	),
    	MinFound = cms.int32( 3 ),
    	hasSelector = cms.vint32( 0, 0 ),
    	TrackProducers = cms.VInputTag( 'hltIter2Merged','hltIter3PFlowTrackSelectionHighPurity' ),
    	LostHitPenalty = cms.double( 20.0 ),
	trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    	newQuality = cms.string( "confirmed" )
    )

    process.HLTIterativeTrackingIteration3 = cms.Sequence( process.hltIter3ClustersRefRemoval + process.hltIter3MaskedMeasurementTrackerEvent + process.hltIter3PixelLayerTriplets + process.hltIter3PFlowPixelTrackingRegions + process.hltIter3PFlowPixelClusterCheck + process.hltIter3PFlowPixelHitDoublets + process.hltIter3PFlowPixelHitTriplets + process.hltIter3PFlowPixelSeeds + process.hltIter3PFlowCkfTrackCandidates + process.hltIter3PFlowCtfWithMaterialTracks + process.hltIter3PFlowTrackCutClassifier + process.hltIter3PFlowTrackSelectionHighPurity )
    #process.HLTIterativeTrackingIter03 = cms.Sequence( process.HLTIterativeTrackingIteration0 + process.HLTIter0TrackAndTauJet4Iter1Sequence + process.HLTIterativeTrackingIteration1 + process.hltIter1Merged + process.HLTIter1TrackAndTauJets4Iter2Sequence + process.HLTIterativeTrackingIteration2 + process.hltIter2Merged + process.HLTIter2TrackAndTauJets4Iter3Sequence )
    process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 + process.HLTIter0TrackAndTauJet4Iter1Sequence + process.HLTIterativeTrackingIteration1 + process.hltIter1Merged + process.HLTIter1TrackAndTauJets4Iter2Sequence + process.HLTIterativeTrackingIteration2 + process.hltIter2Merged + process.HLTIterativeTrackingIteration3 + process.hltIter3Merged  )

    return process
	
def customizeForRecoveryDoublets(process):

    process.hltIter3TrackRefsForJets4Iter4 = cms.EDProducer( "ChargedRefCandidateProducer",
    	src = cms.InputTag( "hltIter3Merged" ),
    	particleType = cms.string( "pi+" )
    )


    process.hltAK4Iter3TrackJets4Iter4 = cms.EDProducer( "FastjetJetProducer",
    	Active_Area_Repeats = cms.int32( 5 ),
    	doAreaFastjet = cms.bool( False ),
    	voronoiRfact = cms.double( 0.9 ),
    	maxBadHcalCells = cms.uint32( 9999999 ),
    	doAreaDiskApprox = cms.bool( False ),
    	maxRecoveredEcalCells = cms.uint32( 9999999 ),
    	jetType = cms.string( "TrackJet" ),
    	minSeed = cms.uint32( 14327 ),
    	Ghost_EtaMax = cms.double( 6.0 ),
    	doRhoFastjet = cms.bool( False ),
    	jetAlgorithm = cms.string( "AntiKt" ),
    	nSigmaPU = cms.double( 1.0 ),
    	GhostArea = cms.double( 0.01 ),
    	Rho_EtaMax = cms.double( 4.4 ),
    	maxBadEcalCells = cms.uint32( 9999999 ),
    	useDeterministicSeed = cms.bool( True ),
    	doPVCorrection = cms.bool( False ),
    	maxRecoveredHcalCells = cms.uint32( 9999999 ),
    	rParam = cms.double( 0.4 ),
    	maxProblematicHcalCells = cms.uint32( 9999999 ),
    	doOutputJets = cms.bool( True ),
    	src = cms.InputTag( "hltIter3TrackRefsForJets4Iter4" ),
    	inputEtMin = cms.double( 0.1 ),
    	puPtMin = cms.double( 0.0 ),
    	srcPVs = cms.InputTag( "hltTrimmedPixelVertices" ),
    	jetPtMin = cms.double( 7.5 ),
    	radiusPU = cms.double( 0.4 ),
    	maxProblematicEcalCells = cms.uint32( 9999999 ),
    	doPUOffsetCorr = cms.bool( False ),
    	inputEMin = cms.double( 0.0 ),
    	useMassDropTagger = cms.bool( False ),
    	muMin = cms.double( -1.0 ),
    	subtractorName = cms.string( "" ),
    	muCut = cms.double( -1.0 ),
    	subjetPtMin = cms.double( -1.0 ),
    	useTrimming = cms.bool( False ),
    	muMax = cms.double( -1.0 ),
    	yMin = cms.double( -1.0 ),
    	useFiltering = cms.bool( False ),
    	rFilt = cms.double( -1.0 ),
    	yMax = cms.double( -1.0 ),
    	zcut = cms.double( -1.0 ),
    	MinVtxNdof = cms.int32( 0 ),
    	MaxVtxZ = cms.double( 30.0 ),
    	UseOnlyVertexTracks = cms.bool( False ),
    	dRMin = cms.double( -1.0 ),
    	nFilt = cms.int32( -1 ),
    	usePruning = cms.bool( False ),
    	maxDepth = cms.int32( -1 ),
    	yCut = cms.double( -1.0 ),
    	DzTrVtxMax = cms.double( 0.5 ),
    	UseOnlyOnePV = cms.bool( True ),
   	rcut_factor = cms.double( -1.0 ),
    	sumRecHits = cms.bool( False ),
    	trimPtFracMin = cms.double( -1.0 ),
    	dRMax = cms.double( -1.0 ),
    	DxyTrVtxMax = cms.double( 0.2 ),
    	useCMSBoostedTauSeedingAlgorithm = cms.bool( False )
    )


    process.hltIter3TrackAndTauJets4Iter4 = cms.EDProducer( "TauJetSelectorForHLTTrackSeeding",
    	fractionMinCaloInTauCone = cms.double( 0.7 ),
    	fractionMaxChargedPUInCaloCone = cms.double( 0.3 ),
    	tauConeSize = cms.double( 0.2 ),
    	ptTrkMaxInCaloCone = cms.double( 1.4 ),
    	isolationConeSize = cms.double( 0.5 ),
    	inputTrackJetTag = cms.InputTag( "hltAK4Iter3TrackJets4Iter4" ),
    	nTrkMaxInCaloCone = cms.int32( 0 ),
    	inputCaloJetTag = cms.InputTag( "hltAK4CaloJetsPFEt5" ),
    	etaMinCaloJet = cms.double( -2.7 ),
    	etaMaxCaloJet = cms.double( 2.7 ),
    	ptMinCaloJet = cms.double( 5.0 ),
    	inputTrackTag = cms.InputTag( "hltIter3Merged" )
    )

    process.hltIter4ClustersRefRemoval = cms.EDProducer( "TrackClusterRemover",
    	trackClassifier = cms.InputTag( '','QualityMasks' ),
    	minNumberOfLayersWithMeasBeforeFiltering = cms.int32( 0 ),
    	maxChi2 = cms.double( 16.0 ),
    	trajectories = cms.InputTag( "hltIter3PFlowTrackSelectionHighPurity" ),
    	oldClusterRemovalInfo = cms.InputTag( "hltIter3ClustersRefRemoval" ),
    	stripClusters = cms.InputTag( "hltSiStripRawToClustersFacility" ),
    	overrideTrkQuals = cms.InputTag( "" ),
    	pixelClusters = cms.InputTag( "hltSiPixelClusters" ),
    	TrackQuality = cms.string( "highPurity" )
    )
    process.hltIter4MaskedMeasurementTrackerEvent = cms.EDProducer( "MaskedMeasurementTrackerEventProducer",
    	clustersToSkip = cms.InputTag( "hltIter4ClustersRefRemoval" ),
    	OnDemand = cms.bool( False ),
    	src = cms.InputTag( "hltSiStripClusters" )
    )

    process.hltIter4PixelLayerPairs = cms.EDProducer( "SeedingLayersEDProducer",
        layerList = cms.vstring( 
		   'BPix1+BPix2', 
		   'BPix1+BPix3', 
		   'BPix2+BPix3',
    		   'BPix1+FPix1_pos',
		   'BPix1+FPix1_neg',
    	           'BPix2+FPix1_pos', 
		   'BPix2+FPix1_neg',
    	           'BPix3+BPix4',
    		   'BPix3+FPix1_pos',
    	           'FPix1_pos+FPix2_pos',
                   'BPix1+BPix4',  
        ),
        MTOB = cms.PSet( ),
        TEC = cms.PSet( ),
        MTID = cms.PSet( ),
        FPix = cms.PSet(
            HitProducer = cms.string( "hltSiPixelRecHits" ),
            hitErrorRZ = cms.double( 0.0036 ),
            useErrorsFromParam = cms.bool( True ),
            TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
            skipClusters = cms.InputTag( "hltIter4ClustersRefRemoval" ),
            hitErrorRPhi = cms.double( 0.0051 )
        ),
        MTEC = cms.PSet( ),
        MTIB = cms.PSet( ),
        TID = cms.PSet( ),
        TOB = cms.PSet( ),
        BPix = cms.PSet(
            HitProducer = cms.string( "hltSiPixelRecHits" ),
            hitErrorRZ = cms.double( 0.006 ),
            useErrorsFromParam = cms.bool( True ),
            TTRHBuilder = cms.string( "hltESPTTRHBuilderPixelOnly" ),
            skipClusters = cms.InputTag( "hltIter4ClustersRefRemoval" ),
            hitErrorRPhi = cms.double( 0.0027 )
        ),
        TIB = cms.PSet( )
    )
    #process.hltIter3PFlowPixelTrackingRegions = cms.EDProducer( "GlobalTrackingRegionFromBeamSpotEDProducer",
    #	RegionPSet = cms.PSet(
    #  		nSigmaZ = cms.double( 0.0 ),
    #  		beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    #  		ptMin = cms.double( 0.9 ),
    #  		originHalfLength = cms.double( 24.0 ),
    #  		originRadius = cms.double( 0.2 ),
    #  		precise = cms.bool( True ),
    #  		useMultipleScattering = cms.bool( False )
   # 	)
    #)


    process.hltIter4PFlowPixelTrackingRegions = cms.EDProducer( "PointSeededTrackingRegionsEDProducer",
    	RegionPSet = cms.PSet( 
     		vertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
    		zErrorVetex = cms.double( 0.1 ),
      		beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
      		zErrorBeamSpot = cms.double( 15.0 ),
      		maxNVertices = cms.int32( 10 ),
      		maxNRegions = cms.int32( 100 ),
      		nSigmaZVertex = cms.double( 3.0 ),
      		nSigmaZBeamSpot = cms.double( 3.0 ),
      		ptMin = cms.double( 0.8 ),
      		mode = cms.string( "VerticesFixed" ),
      		searchOpt = cms.bool( False ),
      		whereToUseMeasurementTracker = cms.string( "ForSiStrips" ),
      		originRadius = cms.double( 0.02 ),
      		measurementTrackerName = cms.InputTag( "hltIter4MaskedMeasurementTrackerEvent" ),
      		precise = cms.bool( True ),
      		deltaEta = cms.double( 1.2 ),
      		deltaPhi = cms.double( 0.5 ),
		points = cms.VPSet(
    		    	cms.PSet(
				eta = cms.double(0.0),
				phi = cms.double(3.0)
	
			)

		)
    	)

    )




    #process.hltIter4PFlowPixelTrackingRegions = cms.EDProducer( "CandidateSeededTrackingRegionsEDProducer",
    #	RegionPSet = cms.PSet( 
     # 	vertexCollection = cms.InputTag( "hltTrimmedPixelVertices" ),
    # 	zErrorVetex = cms.double( 0.05 ),
    #  	beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    #  	zErrorBeamSpot = cms.double( 15.0 ),
    #  	maxNVertices = cms.int32( 10 ),
    #  	maxNRegions = cms.int32( 100 ),
    #  	nSigmaZVertex = cms.double( 3.0 ),
    #  	nSigmaZBeamSpot = cms.double( 3.0 ),
    #  	ptMin = cms.double( 1.2 ),
    #  	mode = cms.string( "VerticesFixed" ),
    #  	input = cms.InputTag( "hltIter3TrackAndTauJets4Iter4" ),
    #  	searchOpt = cms.bool( True ),
    #  	whereToUseMeasurementTracker = cms.string( "ForSiStrips" ),
    #  	originRadius = cms.double( 0.025 ),
    #  	measurementTrackerName = cms.InputTag( "hltIter3MaskedMeasurementTrackerEvent" ),
    #  	precise = cms.bool( True ),
    #  	deltaEta = cms.double( 0.8 ),
    # 	deltaPhi = cms.double( 0.8 )
   # 	)
   # )
    process.hltIter4PFlowPixelClusterCheck = cms.EDProducer( "ClusterCheckerEDProducer",
    	cut = cms.string( "" ),
    	silentClusterCheck = cms.untracked.bool( False ),
    	MaxNumberOfCosmicClusters = cms.uint32( 50000 ),
    	PixelClusterCollectionLabel = cms.InputTag( "hltSiPixelClusters" ),
    	doClusterCheck = cms.bool( False ),
    	MaxNumberOfPixelClusters = cms.uint32( 40000 ),
    	ClusterCollectionLabel = cms.InputTag( "hltSiStripClusters" )
    )
    process.hltIter4PFlowPixelHitDoublets = cms.EDProducer( "HitPairEDProducer",
    	trackingRegions = cms.InputTag( "hltIter4PFlowPixelTrackingRegions" ),
    	layerPairs = cms.vuint32( 0 ),
    	clusterCheck = cms.InputTag( "hltIter4PFlowPixelClusterCheck" ),
    	produceSeedingHitSets = cms.bool( True ),
    	produceIntermediateHitDoublets = cms.bool( False ),
    	maxElement = cms.uint32( 0 ),
    	seedingLayers = cms.InputTag( "hltIter4PixelLayerPairs" )
    )
    process.hltIter4PFlowPixelSeeds = cms.EDProducer( "SeedCreatorFromRegionConsecutiveHitsEDProducer",
    	SeedComparitorPSet = cms.PSet(  ComponentName = cms.string( "none" ) ),
    	forceKinematicWithRegionDirection = cms.bool( False ),
    	magneticField = cms.string( "ParabolicMf" ),
    	SeedMomentumForBOFF = cms.double( 5.0 ),
    	OriginTransverseErrorMultiplier = cms.double( 1.0 ),
    	TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    	MinOneOverPtError = cms.double( 1.0 ),
    	seedingHitSets = cms.InputTag( "hltIter4PFlowPixelHitDoublets" ),
    	propagator = cms.string( "PropagatorWithMaterialParabolicMf" )
    )

    process.HLTIter4GroupedCkfTrajectoryBuilderIT = cms.PSet( 
  	keepOriginalIfRebuildFails = cms.bool( False ),
  	lockHits = cms.bool( True ),
  	maxDPhiForLooperReconstruction = cms.double( 2.0 ),
  	propagatorOpposite = cms.string( "PropagatorWithMaterialParabolicMfOpposite" ),
  	trajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter4PSetTrajectoryFilterIT" ) ),
  	doSeedingRegionRebuilding = cms.bool( False ),
  	useHitsSplitting = cms.bool( False ),
  	maxCand = cms.int32( 2 ),
  	estimator = cms.string( "hltESPChi2ChargeMeasurementEstimator16" ),
  	intermediateCleaning = cms.bool( True ),
  	bestHitOnly = cms.bool( True ),
  	useSameTrajFilter = cms.bool( True ),
  	MeasurementTrackerName = cms.string( "hltESPMeasurementTracker" ),
  	ComponentType = cms.string( "GroupedCkfTrajectoryBuilder" ),
  	lostHitPenalty = cms.double( 30.0 ),
  	requireSeedHitsInRebuild = cms.bool( True ),
  	TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
  	maxPtForLooperReconstruction = cms.double( 0.7 ),
  	cleanTrajectoryAfterInOut = cms.bool( False ),
  	propagatorAlong = cms.string( "PropagatorWithMaterialParabolicMf" ),
  	minNrOfHitsForRebuild = cms.int32( 5 ),
  	alwaysUseInvalidHits = cms.bool( False ),
  	inOutTrajectoryFilter = cms.PSet(  refToPSet_ = cms.string( "HLTIter4PSetTrajectoryFilterIT" ) ),
  	foundHitBonus = cms.double( 5.0 ),
  	updator = cms.string( "hltESPKFUpdator" )
    )

    process.HLTIter4PSetTrajectoryFilterIT = cms.PSet( 
  	minPt = cms.double( 0.3 ),
  	minHitsMinPt = cms.int32( 3 ),
  	ComponentType = cms.string( "CkfBaseTrajectoryFilter" ),
  	maxLostHits = cms.int32( 1 ),
  	maxNumberOfHits = cms.int32( 100 ),
  	maxConsecLostHits = cms.int32( 1 ),
  	minimumNumberOfHits = cms.int32( 4 ),
  	nSigmaMinPt = cms.double( 5.0 ),
  	chargeSignificance = cms.double( -1.0 ),
  	minGoodStripCharge = cms.PSet(  refToPSet_ = cms.string( "HLTSiStripClusterChargeCutNone" ) ),
 	maxCCCLostHits = cms.int32( 0 ),
  	seedExtension = cms.int32( 1 ),
  	strictSeedExtension = cms.bool( False ),
  	minNumberOfHitsForLoopers = cms.int32( 13 ),
  	minNumberOfHitsPerLoop = cms.int32( 4 ),
  	extraNumberOfHitsBeforeTheFirstLoop = cms.int32( 4 ),
  	maxLostHitsFraction = cms.double( 999.0 ),
  	constantValueForLostHitsFractionFilter = cms.double( 1.0 ),
  	seedPairPenalty = cms.int32( 0 ),
  	pixelSeedExtension = cms.bool( False )
    )

    process.hltIter4PFlowCkfTrackCandidates = cms.EDProducer( "CkfTrackCandidateMaker",
	src = cms.InputTag( "hltIter4PFlowPixelSeeds" ),
    	maxSeedsBeforeCleaning = cms.uint32( 1000 ),
    	SimpleMagneticField = cms.string( "ParabolicMf" ),
    	TransientInitialStateEstimatorParameters = cms.PSet( 
     		propagatorAlongTISE = cms.string( "PropagatorWithMaterialParabolicMf" ),
      		numberMeasurementsForFit = cms.int32( 4 ),
      		propagatorOppositeTISE = cms.string( "PropagatorWithMaterialParabolicMfOpposite" )
    	),
    	TrajectoryCleaner = cms.string( "hltESPTrajectoryCleanerBySharedHits" ),
    	MeasurementTrackerEvent = cms.InputTag( "hltIter4MaskedMeasurementTrackerEvent" ),
        cleanTrajectoryAfterInOut = cms.bool( False ),
        useHitsSplitting = cms.bool( False ),
        RedundantSeedCleaner = cms.string( "CachingSeedCleanerBySharedInput" ),
        doSeedingRegionRebuilding = cms.bool( False ),
        maxNSeeds = cms.uint32( 100000 ),
	produceSeedStopReasons = cms.bool( False ),
        TrajectoryBuilderPSet = cms.PSet(  refToPSet_ = cms.string( "HLTIter4PSetTrajectoryBuilderIT" ) ),
    	NavigationSchool = cms.string( "SimpleNavigationSchool" ),
    	TrajectoryBuilder = cms.string( "" )
    )
    process.hltIter4PFlowCtfWithMaterialTracks = cms.EDProducer( "TrackProducer",
    	src = cms.InputTag( "hltIter4PFlowCkfTrackCandidates" ),
    	SimpleMagneticField = cms.string( "ParabolicMf" ),
    	clusterRemovalInfo = cms.InputTag( "" ),
    	beamSpot = cms.InputTag( "hltOnlineBeamSpot" ),
    	MeasurementTrackerEvent = cms.InputTag( "hltIter4MaskedMeasurementTrackerEvent" ),
    	Fitter = cms.string( "hltESPFittingSmootherIT" ),
    	useHitsSplitting = cms.bool( False ),
    	MeasurementTracker = cms.string( "" ),
    	AlgorithmName = cms.string( "hltIter4" ),
    	alias = cms.untracked.string( "ctfWithMaterialTracks" ),
    	NavigationSchool = cms.string( "" ),
    	TrajectoryInEvent = cms.bool( False ),
    	TTRHBuilder = cms.string( "hltESPTTRHBWithTrackAngle" ),
    	GeometricInnerState = cms.bool( True ),
    	useSimpleMF = cms.bool( True ),
    	Propagator = cms.string( "hltESPRungeKuttaTrackerPropagator" )
    )
    process.hltIter4PFlowTrackCutClassifier = cms.EDProducer( "TrackCutClassifier",
    	src = cms.InputTag( "hltIter4PFlowCtfWithMaterialTracks" ),
    	GBRForestLabel = cms.string( "" ),
    	beamspot = cms.InputTag( "hltOnlineBeamSpot" ),
    	vertices = cms.InputTag( "hltTrimmedPixelVertices" ),
    	qualityCuts = cms.vdouble( -0.7, 0.1, 0.7 ),
    	mva = cms.PSet( 
      		minPixelHits = cms.vint32( 0, 0, 0 ),
      		maxDzWrtBS = cms.vdouble( 3.40282346639E38, 24.0, 15.0 ),
      		dr_par = cms.PSet( 
        		d0err = cms.vdouble( 0.003, 0.003, 0.003 ),
        		dr_par2 = cms.vdouble( 3.40282346639E38, 0.3, 0.3 ),
        		dr_par1 = cms.vdouble( 3.40282346639E38, 0.4, 0.4 ),
        		dr_exp = cms.vint32( 4, 4, 4 ),
 		       d0err_par = cms.vdouble( 0.001, 0.001, 0.001 )
		),
      		maxLostLayers = cms.vint32( 1, 1, 1 ),
      		min3DLayers = cms.vint32( 0, 0, 0 ),
      		dz_par = cms.PSet( 
        	dz_par1 = cms.vdouble( 3.40282346639E38, 0.4, 0.4 ),
        	dz_par2 = cms.vdouble( 3.40282346639E38, 0.35, 0.35 ),
        	dz_exp = cms.vint32( 4, 4, 4 )
      	),
      	minNVtxTrk = cms.int32( 3 ),
      	maxDz = cms.vdouble( 0.5, 0.2, 3.40282346639E38 ),
      	minNdof = cms.vdouble( 1.0E-5, 1.0E-5, 1.0E-5 ),
      	maxChi2 = cms.vdouble( 9999.0, 25.0, 16.0 ),
      	maxChi2n = cms.vdouble( 1.2, 1.0, 0.7 ),
      	maxDr = cms.vdouble( 0.5, 0.03, 3.40282346639E38 ),
      	minLayers = cms.vint32( 3, 3, 3 )
    	),
    	GBRForestFileName = cms.string( "" )
    )  
    process.hltIter4PFlowTrackSelectionHighPurity = cms.EDProducer( "TrackCollectionFilterCloner",
    	minQuality = cms.string( "highPurity" ),
    	copyExtras = cms.untracked.bool( True ),
    	copyTrajectories = cms.untracked.bool( False ),
    	originalSource = cms.InputTag( "hltIter4PFlowCtfWithMaterialTracks" ),
    	originalQualVals = cms.InputTag( 'hltIter4PFlowTrackCutClassifier','QualityMasks' ),
    	originalMVAVals = cms.InputTag( 'hltIter4PFlowTrackCutClassifier','MVAValues' )
    )
    process.hltIter4Merged = cms.EDProducer( "TrackListMerger",
    	ShareFrac = cms.double( 0.19 ),
    	writeOnlyTrkQuals = cms.bool( False ),
    	MinPT = cms.double( 0.05 ),
    	allowFirstHitShare = cms.bool( True ),
    	copyExtras = cms.untracked.bool( True ),
    	Epsilon = cms.double( -0.001 ),
    	selectedTrackQuals = cms.VInputTag( 'hltIter3Merged','hltIter4PFlowTrackSelectionHighPurity' ),
    	indivShareFrac = cms.vdouble( 1.0, 1.0 ),
    	MaxNormalizedChisq = cms.double( 1000.0 ),
    	copyMVA = cms.bool( False ),
    	FoundHitBonus = cms.double( 5.0 ),
    	setsToMerge = cms.VPSet( 
      		cms.PSet(  pQual = cms.bool( False ),
        		tLists = cms.vint32( 0, 1 )
      		)
    	),
    	MinFound = cms.int32( 3 ),
    	hasSelector = cms.vint32( 0, 0 ),
    	TrackProducers = cms.VInputTag( 'hltIter3Merged','hltIter4PFlowTrackSelectionHighPurity' ),
    	LostHitPenalty = cms.double( 20.0 ),
    	trackAlgoPriorityOrder = cms.string( "hltESPTrackAlgoPriorityOrder" ),
    	newQuality = cms.string( "confirmed" )
    )
    process.HLTIter3TrackAndTauJets4Iter4Sequence = cms.Sequence( process.hltIter3TrackRefsForJets4Iter4 + process.hltAK4Iter3TrackJets4Iter4 + process.hltIter3TrackAndTauJets4Iter4 )
    process.HLTIterativeTrackingIteration4 = cms.Sequence( process.hltIter4ClustersRefRemoval + process.hltIter4MaskedMeasurementTrackerEvent + process.hltIter4PixelLayerPairs + process.hltIter4PFlowPixelTrackingRegions + process.hltIter4PFlowPixelClusterCheck + process.hltIter4PFlowPixelHitDoublets + process.hltIter4PFlowPixelSeeds + process.hltIter4PFlowCkfTrackCandidates + process.hltIter4PFlowCtfWithMaterialTracks + process.hltIter4PFlowTrackCutClassifier + process.hltIter4PFlowTrackSelectionHighPurity )
    #process.HLTIterativeTrackingIter03 = cms.Sequence( process.HLTIterativeTrackingIteration0 + process.HLTIter0TrackAndTauJet4Iter1Sequence + process.HLTIterativeTrackingIteration1 + process.hltIter1Merged + process.HLTIter1TrackAndTauJets4Iter2Sequence + process.HLTIterativeTrackingIteration2 + process.hltIter2Merged + process.HLTIter2TrackAndTauJets4Iter3Sequence )
    process.HLTIterativeTrackingIter02 = cms.Sequence( process.HLTIterativeTrackingIteration0 + process.HLTIter0TrackAndTauJet4Iter1Sequence + process.HLTIterativeTrackingIteration1 + process.hltIter1Merged + process.HLTIter1TrackAndTauJets4Iter2Sequence + process.HLTIterativeTrackingIteration2 + process.hltIter2Merged + process.HLTIterativeTrackingIteration3 + process.hltIter3Merged  + process.HLTIter3TrackAndTauJets4Iter4Sequence + process.HLTIterativeTrackingIteration4 + process.hltIter4Merged )


    #process.MC_ReducedIterativeTracking_v4 = cms.Path( process.HLTBeginSequence + process.hltPreMCReducedIterativeTracking + process.HLTRecoJetSequenceAK4PrePF + process.HLTDoLocalPixelSequence + process.HLTRecopixelvertexingSequence + process.HLTDoLocalStripSequence + process.HLTIterativeTrackingIter03 + process.HLTEndSequence )
    #process.newSchedule = cms.Schedule(process.MC_ReducedIterativeTracking_v4)
    return process
