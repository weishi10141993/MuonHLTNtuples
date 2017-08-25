
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include "DataFormats/Math/interface/deltaR.h"
#include "MuonHLTNtuples/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

bool           matchMuonWithL3 ( MuonCand, std::vector<HLTMuonCand>);
vector<size_t> matchHltTrackMu ( std::vector<MuonCand>, std::vector<HltTrackCand>, TH1F*, TH1F* );
vector<size_t> matchHltMu ( std::vector<MuonCand>, std::vector<HLTMuonCand> );


void MuonOffHltTrackMatch_IterL3(TString inputfilename="/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/cascade/JSON_08_16/results_alignment/results.root", std::string effmeasured="Cascade"){

  TFile* outfile = TFile::Open(Form("%s_Mu_HltTrack_Match.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  //Create histograms   
  TH1F* DRforOI          = new TH1F("h_DRforOI"          ,"DRforOI"          ,  100,    0,   1. );
  // TH1F* DRforIOL1        = new TH1F("h_DRforIOL1"        ,"DRforIOL1"        ,  100,    0,   1. );
  // TH1F* DRforIOL2        = new TH1F("h_DRforIOL2"        ,"DRforIOL2"        ,  100,    0,   1. );
  TH1F* eta_muons        = new TH1F("h_eta_muon"         , "eta_muon"        , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass   = new TH1F("h_eta_hltOI_pass"   , "eta_hltOI_pass"  , 100,  -2.5,  2.5 );
  // TH1F* eta_hltIOL1_pass = new TH1F("h_eta_hltIOL1_pass" , "eta_hltIOL1_pass", 100,  -2.5,  2.5 );
  // TH1F* eta_hltIOL2_pass  = new TH1F("h_eta_hltIOL2_pass", "eta_hltIOL2_pass", 100,  -2.5,  2.5 );

  TH1F* eta_muons_pass          = new TH1F("h_eta_muon_pass"         , "eta_muon_pass"        , 100,  -2.5,  2.5 );
  TH1F* eta_muons_barrel        = new TH1F("h_eta_muon_barrel"         , "eta_muon_barrel"        , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass_barrel   = new TH1F("h_eta_hltOI_pass_barrel"   , "eta_hltOI_pass_barrel"  , 100,  -2.5,  2.5 );
  TH1F* eta_muons_int           = new TH1F("h_eta_muon_int"            , "eta_muon_int"           , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass_int      = new TH1F("h_eta_hltOI_pass_int"      , "eta_hltOI_pass_int"     , 100,  -2.5,  2.5 );
  TH1F* eta_muons_endcap        = new TH1F("h_eta_muon_endcap"         , "eta_muon_endcap"        , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass_endcap   = new TH1F("h_eta_hltOI_pass_endcap"   , "eta_hltOI_pass_endcap"  , 100,  -2.5,  2.5 );

  TH1F* eta_hltOI   = new TH1F("h_eta_hltOI"   , "eta_hltOI"  , 100,  -2.5,  2.5 );
  TH1F* eta_resta   = new TH1F("h_eta_resta"   , "eta_resta"  , 1000,  -2.5,  2.5 );

  // TH1F* eta_hltIOL1 = new TH1F("h_eta_hltIOL1" , "eta_hltIOL1", 100,  -2.5,  2.5 );
  // TH1F* eta_hltIOL2  = new TH1F("h_eta_hltIOL2", "eta_hltIOL2", 100,  -2.5,  2.5 );

  TH1F* PtHLT         = new TH1F("h_PtHLT"        , "PtHLT"       ,   80,  -0.25,    0.25 );
  TH1F* PtHLT_barrel  = new TH1F("h_PtHLT_barrel" , "PtHLT_barrel",   80,  -0.25,    0.25 );
  TH1F* PtHLT_int     = new TH1F("h_PtHLT_int"    , "PtHLT_int"   ,   80,  -0.25,    0.25 );
  TH1F* PtHLT_endcap  = new TH1F("h_PtHLT_endcap" , "PtHLT_endcap",   80,  -0.25,    0.25 );

  TH1F* PtOI     = new TH1F("h_PtOI"     , "PtOI"     ,   80,  -0.25,    0.25 );
  TH1F* chi2OI   = new TH1F("h_chi2OI"   , "chi2OI"   ,  100,  -10.,   10.);
  TH1F* dxyOI    = new TH1F("h_dxyOI"    , "dxyOI"    ,  100, -0.1,   0.1 );
  TH1F* dzOI     = new TH1F("h_dzOI"     , "dzOI"     ,  100, -0.1,   0.1 );
  TH1F* PixHitOI = new TH1F("h_PixHitOI" , "PixHitOI" ,    25,   -12,     12 );
  TH1F* LayHitOI = new TH1F("h_LayHitOI" , "LayHitOI" ,    21,   -10,     10  );
  TH1F* PixLayOI = new TH1F("h_PixLayOI" , "PixLayOI" ,    21,   -10,     10 );


  //BARREL REGION HISTOGRAMS
  TH1F* PtOI_barrel     = new TH1F("h_PtOI_barrel"     , "PtOI_barrel"     ,   80,  -0.25,    0.25 );
  TH1F* chi2OI_barrel   = new TH1F("h_chi2OI_barrel"   , "chi2OI_barrel"   ,  100,  -10.,   10.);
  TH1F* dxyOI_barrel    = new TH1F("h_dxyOI_barrel"    , "dxyOI_barrel"    ,  100, -0.1,   0.1 );
  TH1F* dzOI_barrel     = new TH1F("h_dzOI_barrel"     , "dzOI_barrel"     ,  100, -0.1,   0.1 );
  TH1F* PixHitOI_barrel = new TH1F("h_PixHitOI_barrel" , "PixHitOI_barrel" ,   25,   -12,     12 );
  TH1F* LayHitOI_barrel = new TH1F("h_LayHitOI_barrel" , "LayHitOI_barrel" ,    21,   -10,     10 );
  TH1F* PixLayOI_barrel = new TH1F("h_PixLayOI_barrel" , "PixLayOI_barrel" ,    21,   -10,     10  );

  // TH1F* PtIOL1_barrel     = new TH1F("h_PtIOL1_barrel"     , "PtIOL1_barrel"     ,   80,  -10,    10 );
  // TH1F* chi2IOL1_barrel   = new TH1F("h_chi2IOL1_barrel"   , "chi2IOL1_barrel"   ,  100,  -10.,   10.);
  // TH1F* dxyIOL1_barrel    = new TH1F("h_dxyIOL1_barrel"    , "dxyIOL1_barrel"    ,  100, -0.1,   0.1 );
  // TH1F* dzIOL1_barrel     = new TH1F("h_dxyIOL1_barrel"    , "dxyIOL1_barrel"    ,  100, -0.1,   0.1 );
  // TH1F* PixHitIOL1_barrel = new TH1F("h_PixHitIOL1_barrel" , "PixHitIOL1_barrel" ,    4,   -2,     2 );
  // TH1F* LayHitIOL1_barrel = new TH1F("h_LayHitIOL1_barrel" , "LayHitIOL1_barrel" ,    4,   -2,     2 );
  // TH1F* PixLayIOL1_barrel = new TH1F("h_PixLayIOL1_barrel" , "PixLayIOL1_barrel" ,    4,   -2,     2 );

  // TH1F* PtIOL2_barrel     = new TH1F("h_PtIOL2_barrel"    , "PtIOL2_barrel"     ,   80,  -10,    10 );
  // TH1F* chi2IOL2_barrel   = new TH1F("h_chi2IOL2_barrel"  , "chi2IOL2_barrel"   ,  100,  -10.,   10.);
  // TH1F* dxyIOL2_barrel    = new TH1F("h_dxyIOL2_barrel"   , "dxyIOL2_barrel"    ,  100, -0.1,   0.1 );
  // TH1F* dzIOL2_barrel     = new TH1F("h_dxyIOL2_barrel"   , "dxyIOL2_barrel"    ,  100, -0.1,   0.1 );
  // TH1F* PixHitIOL2_barrel = new TH1F("h_PixHitIOL2_barrel", "PixHitIOL2_barrel" ,    4,   -2,     2 );
  // TH1F* LayHitIOL2_barrel = new TH1F("h_LayHitIOL2_barrel", "LayHitIOL2_barrel" ,    4,   -2,     2 );
  // TH1F* PixLayIOL2_barrel = new TH1F("h_PixLayIOL2_barrel", "PixLayIOL2_barrel" ,    4,   -2,     2 );

  //Intermediate
  TH1F* PtOI_int     = new TH1F("h_PtOI_int"     , "PtOI_int"     ,   80,  -0.25,    0.25 );
  TH1F* chi2OI_int   = new TH1F("h_chi2OI_int"   , "chi2OI_int"   ,  100,  -10.,   10.);
  TH1F* dxyOI_int    = new TH1F("h_dxyOI_int"    , "dxyOI_int"    ,  100, -0.1,   0.1 );
  TH1F* dzOI_int     = new TH1F("h_dzOI_int"     , "dzOI_int"     ,  100, -0.1,   0.1 );
  TH1F* PixHitOI_int = new TH1F("h_PixHitOI_int" , "PixHitOI_int" ,    25,   -12,     12 );
  TH1F* LayHitOI_int = new TH1F("h_LayHitOI_int" , "LayHitOI_int" ,    21,   -10,     10 );
  TH1F* PixLayOI_int = new TH1F("h_PixLayOI_int" , "PixLayOI_int" ,    21,   -10,     10 );

  // TH1F* PtIOL1_int     = new TH1F("h_PtIOL1_int"     , "PtIOL1_int"     ,   80,  -10,    10 );
  // TH1F* chi2IOL1_int   = new TH1F("h_chi2IOL1_int"   , "chi2IOL1_int"   ,  100,  -10.,   10.);
  // TH1F* dxyIOL1_int    = new TH1F("h_dxyIOL1_int"    , "dxyIOL1_int"    ,  100, -0.1,   0.1 );
  // TH1F* dzIOL1_int     = new TH1F("h_dxyIOL1_int"    , "dxyIOL1_int"    ,  100, -0.1,   0.1 );
  // TH1F* PixHitIOL1_int = new TH1F("h_PixHitIOL1_int" , "PixHitIOL1_int" ,    4,   -2,     2 );
  // TH1F* LayHitIOL1_int = new TH1F("h_LayHitIOL1_int" , "LayHitIOL1_int" ,    4,   -2,     2 );
  // TH1F* PixLayIOL1_int = new TH1F("h_PixLayIOL1_int" , "PixLayIOL1_int" ,    4,   -2,     2 );

  // TH1F* PtIOL2_int     = new TH1F("h_PtIOL2_int"    , "PtIOL2_int"     ,   80,  -10,    10 );
  // TH1F* chi2IOL2_int   = new TH1F("h_chi2IOL2_int"  , "chi2IOL2_int"   ,  100,  -10.,   10.);
  // TH1F* dxyIOL2_int    = new TH1F("h_dxyIOL2_int"   , "dxyIOL2_int"    ,  100, -0.1,   0.1 );
  // TH1F* dzIOL2_int     = new TH1F("h_dxyIOL2_int"   , "dxyIOL2_int"    ,  100, -0.1,   0.1 );
  // TH1F* PixHitIOL2_int = new TH1F("h_PixHitIOL2_int", "PixHitIOL2_int" ,    4,   -2,     2 );
  // TH1F* LayHitIOL2_int = new TH1F("h_LayHitIOL2_int", "LayHitIOL2_int" ,    4,   -2,     2 );
  // TH1F* PixLayIOL2_int = new TH1F("h_PixLayIOL2_int", "PixLayIOL2_int" ,    4,   -2,     2 );


  //endcap
  TH1F* PtOI_endcap     = new TH1F("h_PtOI_endcap"     , "PtOI_endcap"     ,    80,    -0.25,    0.25 );
  TH1F* chi2OI_endcap   = new TH1F("h_chi2OI_endcap"   , "chi2OI_endcap"   ,   100,  -10.,   10.);
  TH1F* dxyOI_endcap    = new TH1F("h_dxyOI_endcap"    , "dxyOI_endcap"    ,   100,  -0.1,   0.1 );
  TH1F* dzOI_endcap     = new TH1F("h_dzOI_endcap"     , "dzOI_endcap"     ,   100,  -0.1,   0.1 );
  TH1F* PixHitOI_endcap = new TH1F("h_PixHitOI_endcap" , "PixHitOI_endcap" ,    25,   -12,     12);
  TH1F* LayHitOI_endcap = new TH1F("h_LayHitOI_endcap" , "LayHitOI_endcap" ,    21,   -10,     10  );
  TH1F* PixLayOI_endcap = new TH1F("h_PixLayOI_endcap" , "PixLayOI_endcap" ,    21,   -10,     10  );

  // TH1F* PtIOL1_endcap     = new TH1F("h_PtIOL1_endcap"     , "PtIOL1_endcap"     ,   80,  -10,    10 );
  // TH1F* chi2IOL1_endcap   = new TH1F("h_chi2IOL1_endcap"   , "chi2IOL1_endcap"   ,  100,  -10.,   10.);
  // TH1F* dxyIOL1_endcap    = new TH1F("h_dxyIOL1_endcap"    , "dxyIOL1_endcap"    ,  100, -0.1,   0.1 );
  // TH1F* dzIOL1_endcap     = new TH1F("h_dxyIOL1_endcap"    , "dxyIOL1_endcap"    ,  100, -0.1,   0.1 );
  // TH1F* PixHitIOL1_endcap = new TH1F("h_PixHitIOL1_endcap" , "PixHitIOL1_endcap" ,    4,   -2,     2 );
  // TH1F* LayHitIOL1_endcap = new TH1F("h_LayHitIOL1_endcap" , "LayHitIOL1_endcap" ,    4,   -2,     2 );
  // TH1F* PixLayIOL1_endcap = new TH1F("h_PixLayIOL1_endcap" , "PixLayIOL1_endcap" ,    4,   -2,     2 );

  // TH1F* PtIOL2_endcap     = new TH1F("h_PtIOL2_endcap"    , "PtIOL2_endcap"     ,   80,  -10,    10 );
  // TH1F* chi2IOL2_endcap   = new TH1F("h_chi2IOL2_endcap"  , "chi2IOL2_endcap"   ,  100,  -10.,   10.);
  // TH1F* dxyIOL2_endcap    = new TH1F("h_dxyIOL2_endcap"   , "dxyIOL2_endcap"    ,  100, -0.1,   0.1 );
  // TH1F* dzIOL2_endcap     = new TH1F("h_dxyIOL2_endcap"   , "dxyIOL2_endcap"    ,  100, -0.1,   0.1 );
  // TH1F* PixHitIOL2_endcap = new TH1F("h_PixHitIOL2_endcap", "PixHitIOL2_endcap" ,    4,   -2,     2 );
  // TH1F* LayHitIOL2_endcap = new TH1F("h_LayHitIOL2_endcap", "LayHitIOL2_endcap" ,    4,   -2,     2 );
  // TH1F* PixLayIOL2_endcap = new TH1F("h_PixLayIOL2_endcap", "PixLayIOL2_endcap" ,    4,   -2,     2 );

  
  TFile* inputfile = TFile::Open(inputfilename, "READ"); 
  std::cout << "input file: " << inputfile -> GetName() << std::endl;

  TTree *tree = (TTree*) inputfile -> Get("muonNtuples/muonTree"); 
  
  if (!tree) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }
    
  MuonEvent* ev      = new MuonEvent(); 
  TBranch*  evBranch = tree->GetBranch("event");
  evBranch -> SetAddress(&ev);

  int nentries = tree->GetEntriesFast();
  std::cout << "Number of entries = " << nentries << std::endl;

  for (Int_t eventNo=0; eventNo < nentries; eventNo++) {

    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    
    unsigned int nmuons   = ev -> muons.size();
    unsigned int nhltOI   = ev -> hltTrackOI.size(); 

    // unsigned int nhltIOL1 = ev -> hltTrackIOL1.size();
    // unsigned int nhltIOL2 = ev -> hltTrackIOL2.size();

    vector<size_t> matchTrackOI = matchHltTrackMu(ev -> muons, ev -> hltTrackOI, DRforOI, eta_hltOI);
    vector<size_t> matchTKmuon  = matchHltMu(ev -> muons, ev -> tkmuons);

    // vector< size_t> matchTrackIOL1 = matchHltTrackMu (ev -> muons, ev -> hltTrackIOL1, DRforIOL1, eta_hltIOL1);
    // vector< size_t> matchTrackIOL2 = matchHltTrackMu (ev -> muons, ev -> hltTrackIOL2, DRforIOL2, eta_hltIOL2);

    for (size_t imu = 0; imu < nmuons; imu++){

      //if ( !matchMuonWithL3(ev->muons.at(imu),ev -> tkmuons) ) continue;

      int ihltOI = matchTrackOI[imu];
      int itkmu  = matchTKmuon[imu];
      // int ihltIOL1 = matchTrackIOL1[imu];
      // int ihltIOL2 = matchTrackIOL2[imu];

      eta_muons -> Fill (ev -> muons.at(imu).eta);

      if ( itkmu > -1){

	PtHLT -> Fill (((ev -> muons.at(imu).innerpt) - (ev -> tkmuons.at(itkmu).trkpt))/(ev -> muons.at(imu).innerpt));

	if (fabs(ev -> muons.at(imu).eta) <= 0.9 ){
	  PtHLT_barrel -> Fill (((ev -> muons.at(imu).innerpt) - (ev -> tkmuons.at(itkmu).trkpt))/(ev -> muons.at(imu).innerpt));
	}

	if (fabs(ev -> muons.at(imu).eta)>0.9 && fabs(ev -> muons.at(imu).eta)<1.6){
	  PtHLT_int -> Fill (((ev -> muons.at(imu).innerpt) - (ev -> tkmuons.at(itkmu).trkpt))/(ev -> muons.at(imu).innerpt));
	}

	if ( fabs(ev -> muons.at(imu).eta)>=1.6){
	  PtHLT_endcap -> Fill (((ev -> muons.at(imu).innerpt) - (ev -> tkmuons.at(itkmu).trkpt))/(ev -> muons.at(imu).innerpt));
	}
      }	

      if ( ihltOI > -1){

	eta_muons_pass -> Fill (ev -> muons.at(imu).eta);
	eta_hltOI_pass -> Fill (ev -> hltTrackOI.at(ihltOI).eta);
	eta_resta      -> Fill (ev -> muons.at(imu).eta - ev -> hltTrackOI.at(ihltOI).eta);
	PtOI     -> Fill (((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackOI.at(ihltOI).pt))/(ev -> muons.at(imu).innerpt));
	chi2OI   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackOI.at(ihltOI).chi2));
	dxyOI    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackOI.at(ihltOI).dxy));
	dzOI     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackOI.at(ihltOI).dz));
	PixHitOI -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackOI.at(ihltOI).pixelHits));
	LayHitOI -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackOI.at(ihltOI).layerHits));
	PixLayOI -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackOI.at(ihltOI).pixelLayers));


	if (fabs(ev -> muons.at(imu).eta) <= 0.9 ){
	  eta_muons_barrel -> Fill (ev -> muons.at(imu).eta);
	  eta_hltOI_pass_barrel -> Fill (ev -> hltTrackOI.at(ihltOI).eta);
	  PtOI_barrel     -> Fill (((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackOI.at(ihltOI).pt))/(ev -> muons.at(imu).innerpt));
	  chi2OI_barrel   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackOI.at(ihltOI).chi2));
	  dxyOI_barrel    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackOI.at(ihltOI).dxy));
	  dzOI_barrel     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackOI.at(ihltOI).dz));
	  PixHitOI_barrel -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackOI.at(ihltOI).pixelHits));
	  LayHitOI_barrel -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackOI.at(ihltOI).layerHits));
	  PixLayOI_barrel -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackOI.at(ihltOI).pixelLayers));
	} 

	// if ( ihltIOL1 > -1){
	//   eta_hltIOL1_pass -> Fill (ev -> hltTrackIOL1.at(ihltIOL1).eta);
	//   PtIOL1_barrel     -> Fill ((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackIOL1.at(ihltIOL1).pt));
	//   chi2IOL1_barrel   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackIOL1.at(ihltIOL1).chi2));
	//   dxyIOL1_barrel    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackIOL1.at(ihltIOL1).dxy));
	//   dzIOL1_barrel     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackIOL1.at(ihltIOL1).dz));
	//   PixHitIOL1_barrel -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackIOL1.at(ihltIOL1).pixelHits));
	//   LayHitIOL1_barrel -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackIOL1.at(ihltIOL1).layerHits));
	//   PixLayIOL1_barrel -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackIOL1.at(ihltIOL1).pixelLayers));
	// } 

	// if ( ihltIOL2 > -1){
	//   //std::cout << " pasa el L2" << endl;
	//   eta_hltIOL2_pass -> Fill (ev -> hltTrackIOL2.at(ihltIOL2).eta);
	//   PtIOL2_barrel     -> Fill ((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackIOL2.at(ihltIOL2).pt));
	//   chi2IOL2_barrel   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackIOL2.at(ihltIOL2).chi2));
	//   dxyIOL2_barrel    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackIOL2.at(ihltIOL2).dxy));
	//   dzIOL2_barrel     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackIOL2.at(ihltIOL2).dz));
	//   PixHitIOL2_barrel -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackIOL2.at(ihltIOL2).pixelHits));
	//   LayHitIOL2_barrel -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackIOL2.at(ihltIOL2).layerHits));
	//   PixLayIOL2_barrel -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackIOL2.at(ihltIOL2).pixelLayers));
	// } 
      
      
	if (fabs(ev -> muons.at(imu).eta)>0.9 && fabs(ev -> muons.at(imu).eta)<1.6){

	  eta_muons_int -> Fill (ev -> muons.at(imu).eta);
	  eta_hltOI_pass_int -> Fill (ev -> hltTrackOI.at(ihltOI).eta);

	  PtOI_int     -> Fill (((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackOI.at(ihltOI).pt))/((ev -> muons.at(imu).innerpt)));
	  chi2OI_int   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackOI.at(ihltOI).chi2));
	  dxyOI_int    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackOI.at(ihltOI).dxy));
	  dzOI_int     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackOI.at(ihltOI).dz));
	  PixHitOI_int -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackOI.at(ihltOI).pixelHits));
	  LayHitOI_int -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackOI.at(ihltOI).layerHits));
	  PixLayOI_int -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackOI.at(ihltOI).pixelLayers));
	} 

	// if ( ihltIOL1 > -1){
	//   eta_hltIOL1_pass -> Fill (ev -> hltTrackIOL1.at(ihltIOL1).eta);
	//   PtIOL1_int     -> Fill ((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackIOL1.at(ihltIOL1).pt));
	//   chi2IOL1_int   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackIOL1.at(ihltIOL1).chi2));
	//   dxyIOL1_int    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackIOL1.at(ihltIOL1).dxy));
	//   dzIOL1_int     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackIOL1.at(ihltIOL1).dz));
	//   PixHitIOL1_int -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackIOL1.at(ihltIOL1).pixelHits));
	//   LayHitIOL1_int -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackIOL1.at(ihltIOL1).layerHits));
	//   PixLayIOL1_int -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackIOL1.at(ihltIOL1).pixelLayers));
	// } 

	// if ( ihltIOL2 > -1){
	//   eta_hltIOL2_pass -> Fill (ev -> hltTrackIOL2.at(ihltIOL2).eta);
	//   PtIOL2_int     -> Fill ((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackIOL2.at(ihltIOL2).pt));
	//   chi2IOL2_int   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackIOL2.at(ihltIOL2).chi2));
	//   dxyIOL2_int    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackIOL2.at(ihltIOL2).dxy));
	//   dzIOL2_int     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackIOL2.at(ihltIOL2).dz));
	//   PixHitIOL2_int -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackIOL2.at(ihltIOL2).pixelHits));
	//   LayHitIOL2_int -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackIOL2.at(ihltIOL2).layerHits));
	//   PixLayIOL2_int -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackIOL2.at(ihltIOL2).pixelLayers));
	// }

	if ( fabs(ev -> muons.at(imu).eta) >= 1.6){
	  eta_muons_endcap -> Fill (ev -> muons.at(imu).eta);
	  eta_hltOI_pass_endcap -> Fill (ev -> hltTrackOI.at(ihltOI).eta);
	  PtOI_endcap     -> Fill (((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackOI.at(ihltOI).pt))/(ev -> muons.at(imu).innerpt));
	  chi2OI_endcap   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackOI.at(ihltOI).chi2));
	  dxyOI_endcap    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackOI.at(ihltOI).dxy));
	  dzOI_endcap     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackOI.at(ihltOI).dz));
	  PixHitOI_endcap -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackOI.at(ihltOI).pixelHits));
	  LayHitOI_endcap -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackOI.at(ihltOI).layerHits));
	  PixLayOI_endcap -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackOI.at(ihltOI).pixelLayers));
	}

	// if ( ihltIOL1 > -1){
	//   eta_hltIOL1_pass -> Fill (ev -> hltTrackIOL1.at(ihltIOL1).eta);
	//   PtIOL1_endcap     -> Fill ((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackIOL1.at(ihltIOL1).pt));
	//   chi2IOL1_endcap   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackIOL1.at(ihltIOL1).chi2));
	//   dxyIOL1_endcap    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackIOL1.at(ihltIOL1).dxy));
	//   dzIOL1_endcap     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackIOL1.at(ihltIOL1).dz));
	//   PixHitIOL1_endcap -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackIOL1.at(ihltIOL1).pixelHits));
	//   LayHitIOL1_endcap -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackIOL1.at(ihltIOL1).layerHits));
	//   PixLayIOL1_endcap -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackIOL1.at(ihltIOL1).pixelLayers));
	// } 

	// if ( ihltIOL2 > -1){
	//   eta_hltIOL2_pass -> Fill (ev -> hltTrackIOL2.at(ihltIOL2).eta);
	//   PtIOL2_endcap     -> Fill ((ev -> muons.at(imu).innerpt)    - (ev -> hltTrackIOL2.at(ihltIOL2).pt));
	//   chi2IOL2_endcap   -> Fill ((ev -> muons.at(imu).innerchi2)  - (ev -> hltTrackIOL2.at(ihltIOL2).chi2));
	//   dxyIOL2_endcap    -> Fill ((ev -> muons.at(imu).innerdxy)   - (ev -> hltTrackIOL2.at(ihltIOL2).dxy));
	//   dzIOL2_endcap     -> Fill ((ev -> muons.at(imu).innerdz)    - (ev -> hltTrackIOL2.at(ihltIOL2).dz));
	//   PixHitIOL2_endcap -> Fill ((ev -> muons.at(imu).innerpixelHits)  - (ev -> hltTrackIOL2.at(ihltIOL2).pixelHits));
	//   LayHitIOL2_endcap -> Fill ((ev -> muons.at(imu).innerlayerHits)  - (ev -> hltTrackIOL2.at(ihltIOL2).layerHits));
	//   PixLayIOL2_endcap -> Fill ((ev -> muons.at(imu).innerpixelLayers)- (ev -> hltTrackIOL2.at(ihltIOL2).pixelLayers));
	// }

      } 
    }
  }


  
  //Writing the histograms in a file.
  outfile      -> cd();

  DRforOI      -> Write();
  // DRforIOL1    -> Write();
  // DRforIOL2    -> Write();
  eta_hltOI    -> Write();
  // eta_hltIOL1  -> Write();
  // eta_hltIOL2  -> Write();

  eta_muons        -> Write();
  eta_hltOI_pass   -> Write();
  eta_resta        -> Write();


  eta_muons_pass       -> Write();
  eta_muons_barrel     -> Write();
  eta_hltOI_pass_barrel-> Write();
  eta_muons_int        -> Write();
  eta_hltOI_pass_int   -> Write();
  eta_muons_endcap     -> Write();
  eta_hltOI_pass_endcap-> Write();
  // eta_hltIOL1_pass -> Write();
  // eta_hltIOL2_pass -> Write();

  PtHLT        -> Write();
  PtHLT_barrel -> Write();
  PtHLT_int    -> Write();
  PtHLT_endcap -> Write();


  PtOI       -> Write();
  chi2OI     -> Write();
  dxyOI      -> Write();
  dzOI       -> Write();
  PixHitOI   -> Write();
  LayHitOI   -> Write();
  PixLayOI   -> Write();

  //Write histograms barrel region
  PtOI_barrel       -> Write();
  chi2OI_barrel     -> Write();
  dxyOI_barrel      -> Write();
  dzOI_barrel       -> Write();
  PixHitOI_barrel   -> Write();
  LayHitOI_barrel   -> Write();
  PixLayOI_barrel   -> Write();
  
  // PtIOL1_barrel     -> Write();
  // chi2IOL1_barrel   -> Write();
  // dxyIOL1_barrel    -> Write();
  // dzIOL1_barrel     -> Write();
  // PixHitIOL1_barrel -> Write();
  // LayHitIOL1_barrel -> Write();
  // PixLayIOL1_barrel -> Write();

  // PtIOL2_barrel     -> Write();
  // chi2IOL2_barrel   -> Write();
  // dxyIOL2_barrel    -> Write();
  // dzIOL2_barrel     -> Write();
  // PixHitIOL2_barrel -> Write();
  // LayHitIOL2_barrel -> Write();
  // PixLayIOL2_barrel -> Write();
  
  //Write histograms overlap
  PtOI_int      -> Write();
  chi2OI_int    -> Write();
  dxyOI_int     -> Write();
  dzOI_int      -> Write();
  PixHitOI_int  -> Write();
  LayHitOI_int  -> Write();
  PixLayOI_int  -> Write();

  // PtIOL1_int     -> Write();
  // chi2IOL1_int   -> Write();
  // dxyIOL1_int    -> Write();
  // dzIOL1_int     -> Write();
  // PixHitIOL1_int -> Write();
  // LayHitIOL1_int -> Write();
  // PixLayIOL1_int -> Write();

  // PtIOL2_int     -> Write();
  // chi2IOL2_int   -> Write();
  // dxyIOL2_int    -> Write();
  // dzIOL2_int     -> Write();
  // PixHitIOL2_int -> Write();
  // LayHitIOL2_int -> Write();
  // PixLayIOL2_int -> Write();

  //Write histograms endcap region
  PtOI_endcap       -> Write();
  chi2OI_endcap     -> Write();
  dxyOI_endcap      -> Write();
  dzOI_endcap       -> Write();
  PixHitOI_endcap   -> Write();
  LayHitOI_endcap   -> Write();
  PixLayOI_endcap   -> Write();

  // PtIOL1_endcap     -> Write();
  // chi2IOL1_endcap   -> Write();
  // dxyIOL1_endcap    -> Write();
  // dzIOL1_endcap     -> Write();
  // PixHitIOL1_endcap -> Write();
  // LayHitIOL1_endcap -> Write();
  // PixLayIOL1_endcap -> Write();

  // PtIOL2_endcap     -> Write();
  // chi2IOL2_endcap   -> Write();
  // dxyIOL2_endcap    -> Write();
  // dzIOL2_endcap     -> Write();
  // PixHitIOL2_endcap -> Write();
  // LayHitIOL2_endcap -> Write();
  // PixLayIOL2_endcap -> Write();


  outfile      -> Close();  
  
  return;

}


bool matchMuonWithL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

  bool match = false;
  float minDR = 0.1;
  float theDR = 100;
  for ( std::vector<HLTMuonCand>::const_iterator it = L3cands.begin(); it != L3cands.end(); ++it ) { 
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi); 
    if (theDR < minDR){ 
      minDR = theDR;
      match = true;
    }
  }
  return match;
}


vector<size_t> matchHltTrackMu (std::vector<MuonCand> mu, std::vector<HltTrackCand> hltTrack, TH1F* DeltaRplot, TH1F* eta_plot){

  unsigned int nmuons   = mu.size();
  unsigned int nhlt     = hltTrack.size(); 
 
  float deltaRMat[nmuons][nhlt];
  double theDeltaR;
 
  for ( int ihlt=0; ihlt<hltTrack.size(); ihlt++ ) {
    eta_plot -> Fill ( hltTrack.at(ihlt).eta );
  }
  
  //create matrix with all possible deltaR and fill plot
  for (int imu = 0; imu < nmuons; imu++){ 
    for ( int ihlt=0; ihlt < hltTrack.size(); ihlt++ ) {
      deltaRMat[imu][ihlt]=deltaR(mu.at(imu).eta, mu.at(imu).phi, hltTrack.at(ihlt).eta, hltTrack.at(ihlt).phi);
      theDeltaR = deltaRMat[imu][ihlt];
      DeltaRplot -> Fill(theDeltaR);
    }
  }

  float maxDeltaR = 0.1;
  std::vector<size_t>  match(nmuons, -1);


  for (size_t k=0; k < nmuons; k++){
    double minDeltaR = maxDeltaR;
    size_t muon = -999;
    size_t hlt = -999;

    for(int imu = 0; imu < nmuons; imu++){
      for ( int ihlt=0; ihlt < hltTrack.size(); ihlt++ ) {
	if(deltaRMat[imu][ihlt] < minDeltaR){

	  minDeltaR = deltaRMat[imu][ihlt];
	  hlt = ihlt;
	  muon = imu;

	}
      }
    }

    if(minDeltaR < maxDeltaR){
      match[muon] = hlt;
      
      for(int imu = 0;imu < nmuons; imu++){
	for(int ihlt = 0; ihlt < nhlt; ihlt++){
	  if (imu == muon || ihlt == hlt) deltaRMat[imu][ihlt] = 1000;
	}
      }
    } 
  }
  return match;

}



vector<size_t> matchHltMu (std::vector<MuonCand> mu, std::vector<HLTMuonCand> hltmuon){

  unsigned int nmuons   = mu.size();
  unsigned int nhlt     = hltmuon.size(); 
 
  float deltaRMat[nmuons][nhlt];
  double theDeltaR;
  
  //create matrix with all possible deltaR and fill plot
  for (int imu = 0; imu < nmuons; imu++){ 
    for ( int ihlt=0; ihlt < hltmuon.size(); ihlt++ ) {
      deltaRMat[imu][ihlt]=deltaR(mu.at(imu).eta, mu.at(imu).phi, hltmuon.at(ihlt).eta, hltmuon.at(ihlt).phi);
      theDeltaR = deltaRMat[imu][ihlt];
    }
  }

  float maxDeltaR = 0.1;
  std::vector<size_t>  match(nmuons, -1);


  for (size_t k=0; k < nmuons; k++){
    double minDeltaR = maxDeltaR;
    size_t muon = -999;
    size_t hlt = -999;

    for(int imu = 0; imu < nmuons; imu++){
      for ( int ihlt=0; ihlt < hltmuon.size(); ihlt++ ) {
	if(deltaRMat[imu][ihlt] < minDeltaR){

	  minDeltaR = deltaRMat[imu][ihlt];
	  hlt = ihlt;
	  muon = imu;

	}
      }
    }

    if(minDeltaR < maxDeltaR){
      match[muon] = hlt;
      
      for(int imu = 0;imu < nmuons; imu++){
	for(int ihlt = 0; ihlt < nhlt; ihlt++){
	  if (imu == muon || ihlt == hlt) deltaRMat[imu][ihlt] = 1000;
	}
      }
    } 
  }
  return match;

}

