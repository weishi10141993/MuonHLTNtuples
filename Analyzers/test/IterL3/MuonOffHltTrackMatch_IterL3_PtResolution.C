
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
void printProgBar(int);

void MuonOffHltTrackMatch_IterL3_PtResolution(TString inputfilename="/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/cascade/JSON_08_16/results_alignment/results.root", std::string effmeasured="Cascade"){

  TFile* outfile = TFile::Open(Form("%s_Mu_HltTrack_Match.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  //Create histograms   
  TH1F* DRforOI          = new TH1F("h_DRforOI"          ,"DRforOI"          ,  100,    0,   1. );
  TH1F* eta_muons        = new TH1F("h_eta_muon"         , "eta_muon"        , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass   = new TH1F("h_eta_hltOI_pass"   , "eta_hltOI_pass"  , 100,  -2.5,  2.5 );

  TH1F* eta_muons_pass          = new TH1F("h_eta_muon_pass"         , "eta_muon_pass"        , 100,  -2.5,  2.5 );
  TH1F* eta_muons_barrel        = new TH1F("h_eta_muon_barrel"         , "eta_muon_barrel"        , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass_barrel   = new TH1F("h_eta_hltOI_pass_barrel"   , "eta_hltOI_pass_barrel"  , 100,  -2.5,  2.5 );
  TH1F* eta_muons_int           = new TH1F("h_eta_muon_int"            , "eta_muon_int"           , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass_int      = new TH1F("h_eta_hltOI_pass_int"      , "eta_hltOI_pass_int"     , 100,  -2.5,  2.5 );
  TH1F* eta_muons_endcap        = new TH1F("h_eta_muon_endcap"         , "eta_muon_endcap"        , 100,  -2.5,  2.5 );
  TH1F* eta_hltOI_pass_endcap   = new TH1F("h_eta_hltOI_pass_endcap"   , "eta_hltOI_pass_endcap"  , 100,  -2.5,  2.5 );

  TH1F* eta_hltOI   = new TH1F("h_eta_hltOI"   , "eta_hltOI"  , 100,  -2.5,  2.5 );
  TH1F* eta_resta   = new TH1F("h_eta_resta"   , "eta_resta"  , 1000,  -2.5,  2.5 );

  TH1F* PtHLT         = new TH1F("h_PtHLT"        , "PtHLT"       ,   80,  -0.25,    0.25 );
  TH1F* PtHLT_barrel  = new TH1F("h_PtHLT_barrel" , "PtHLT_barrel",   80,  -0.25,    0.25 );
  TH1F* PtHLT_int     = new TH1F("h_PtHLT_int"    , "PtHLT_int"   ,   80,  -0.25,    0.25 );
  TH1F* PtHLT_endcap  = new TH1F("h_PtHLT_endcap" , "PtHLT_endcap",   80,  -0.25,    0.25 );

  TH1F* QoverPtHLT         = new TH1F("h_QoverPtHLT"        , "QoverPtHLT"       ,   80,  -0.25,    0.25 );
  TH1F* QoverPtHLT_barrel  = new TH1F("h_QoverPtHLT_barrel" , "QoverPtHLT_barrel",   80,  -0.25,    0.25 );
  TH1F* QoverPtHLT_int     = new TH1F("h_QoverPtHLT_int"    , "QoverPtHLT_int"   ,   80,  -0.25,    0.25 );
  TH1F* QoverPtHLT_endcap  = new TH1F("h_QoverPtHLT_endcap" , "QoverPtHLT_endcap",   80,  -0.25,    0.25 );

  TH1F* PtHLTGlb         = new TH1F("h_PtHLTGlb"        , "PtHLTGlb"       ,   80,  -0.25,    0.25 );
  TH1F* PtHLTGlb_barrel  = new TH1F("h_PtHLTGlb_barrel" , "PtHLTGlb_barrel",   80,  -0.25,    0.25 );
  TH1F* PtHLTGlb_int     = new TH1F("h_PtHLTGlb_int"    , "PtHLTGlb_int"   ,   80,  -0.25,    0.25 );
  TH1F* PtHLTGlb_endcap  = new TH1F("h_PtHLTGlb_endcap" , "PtHLTGlb_endcap",   80,  -0.25,    0.25 );

  TH1F* QoverPtHLTGlb         = new TH1F("h_QoverPtHLTGlb"        , "QoverPtHLTGlb"       ,   80,  -0.25,    0.25 );
  TH1F* QoverPtHLTGlb_barrel  = new TH1F("h_QoverPtHLTGlb_barrel" , "QoverPtHLTGlb_barrel",   80,  -0.25,    0.25 );
  TH1F* QoverPtHLTGlb_int     = new TH1F("h_QoverPtHLTGlb_int"    , "QoverPtHLTGlb_int"   ,   80,  -0.25,    0.25 );
  TH1F* QoverPtHLTGlb_endcap  = new TH1F("h_QoverPtHLTGlb_endcap" , "QoverPtHLTGlb_endcap",   80,  -0.25,    0.25 );

  TH1F* PtL3OI         = new TH1F("h_PtL3OI"        , "PtL3OI"       ,   80,  -0.25,    0.25 );
  TH1F* PtL3OI_barrel  = new TH1F("h_PtL3OI_barrel" , "PtL3OI_barrel",   80,  -0.25,    0.25 );
  TH1F* PtL3OI_int     = new TH1F("h_PtL3OI_int"    , "PtL3OI_int"   ,   80,  -0.25,    0.25 );
  TH1F* PtL3OI_endcap  = new TH1F("h_PtL3OI_endcap" , "PtL3OI_endcap",   80,  -0.25,    0.25 );

  TH1F* QoverPtL3OI         = new TH1F("h_QoverPtL3OI"        , "QoverPtL3OI"       ,   80,  -0.25,    0.25 );
  TH1F* QoverPtL3OI_barrel  = new TH1F("h_QoverPtL3OI_barrel" , "QoverPtL3OI_barrel",   80,  -0.25,    0.25 );
  TH1F* QoverPtL3OI_int     = new TH1F("h_QoverPtL3OI_int"    , "QoverPtL3OI_int"   ,   80,  -0.25,    0.25 );
  TH1F* QoverPtL3OI_endcap  = new TH1F("h_QoverPtL3OI_endcap" , "QoverPtL3OI_endcap",   80,  -0.25,    0.25 );

  TH1F* PtL3OIGlb         = new TH1F("h_PtL3OIGlb"        , "PtL3OIGlb"       ,   80,  -0.25,    0.25 );
  TH1F* PtL3OIGlb_barrel  = new TH1F("h_PtL3OIGlb_barrel" , "PtL3OIGlb_barrel",   80,  -0.25,    0.25 );
  TH1F* PtL3OIGlb_int     = new TH1F("h_PtL3OIGlb_int"    , "PtL3OIGlb_int"   ,   80,  -0.25,    0.25 );
  TH1F* PtL3OIGlb_endcap  = new TH1F("h_PtL3OIGlb_endcap" , "PtL3OIGlb_endcap",   80,  -0.25,    0.25 );

  TH1F* QoverPtL3OIGlb         = new TH1F("h_QoverPtL3OIGlb"        , "QoverPtL3OIGlb"       ,   80,  -0.25,    0.25 );
  TH1F* QoverPtL3OIGlb_barrel  = new TH1F("h_QoverPtL3OIGlb_barrel" , "QoverPtL3OIGlb_barrel",   80,  -0.25,    0.25 );
  TH1F* QoverPtL3OIGlb_int     = new TH1F("h_QoverPtL3OIGlb_int"    , "QoverPtL3OIGlb_int"   ,   80,  -0.25,    0.25 );
  TH1F* QoverPtL3OIGlb_endcap  = new TH1F("h_QoverPtL3OIGlb_endcap" , "QoverPtL3OIGlb_endcap",   80,  -0.25,    0.25 );
  
  TH1F* PtOI          = new TH1F("h_PtOI"         , "PtOI"        ,   80,  -0.25,    0.25 );
  TH1F* PtOI_barrel   = new TH1F("h_PtOI_barrel"  , "PtOI_barrel" ,   80,  -0.25,    0.25 );
  TH1F* PtOI_int      = new TH1F("h_PtOI_int"     , "PtOI_int"    ,   80,  -0.25,    0.25 );
  TH1F* PtOI_endcap   = new TH1F("h_PtOI_endcap"  , "PtOI_endcap" ,   80,  -0.25,    0.25 );

  TH1F* QoverPtOI          = new TH1F("h_QoverPtOI"         , "QoverPtOI"        ,   80,  -0.25,    0.25 );
  TH1F* QoverPtOI_barrel   = new TH1F("h_QoverPtOI_barrel"  , "QoverPtOI_barrel" ,   80,  -0.25,    0.25 );
  TH1F* QoverPtOI_int      = new TH1F("h_QoverPtOI_int"     , "QoverPtOI_int"    ,   80,  -0.25,    0.25 );
  TH1F* QoverPtOI_endcap   = new TH1F("h_QoverPtOI_endcap"  , "QoverPtOI_endcap" ,   80,  -0.25,    0.25 );
 
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
    printProgBar((int)(eventNo*100./nentries));

    unsigned int nmuons   = ev -> muons.size();
    unsigned int nhltOI   = ev -> hltTrackOI.size(); 

    // unsigned int nhltIOL1 = ev -> hltTrackIOL1.size();
    // unsigned int nhltIOL2 = ev -> hltTrackIOL2.size();

    vector<size_t> matchTrackOI = matchHltTrackMu(ev -> muons, ev -> hltTrackOI, DRforOI, eta_hltOI);
    vector<size_t> matchTKmuon  = matchHltMu(ev -> muons, ev -> tkmuons);
    vector<size_t> matchHLTmuon = matchHltMu(ev -> muons, ev -> hltmuons);

    // vector< size_t> matchTrackIOL1 = matchHltTrackMu (ev -> muons, ev -> hltTrackIOL1, DRforIOL1, eta_hltIOL1);
    // vector< size_t> matchTrackIOL2 = matchHltTrackMu (ev -> muons, ev -> hltTrackIOL2, DRforIOL2, eta_hltIOL2);

    for (size_t imu = 0; imu < nmuons; imu++){

      //if ( !matchMuonWithL3(ev->muons.at(imu),ev -> tkmuons) ) continue;

      int ihltOI = matchTrackOI[imu];
      int itkmu  = matchTKmuon[imu];
      int ihltmu = matchHLTmuon[imu];
      
      float mueta   = ev -> muons.at(imu).eta;
      eta_muons -> Fill (mueta);
      
      float mupttrk = ev -> muons.at(imu).innerpt;
      float muptglb = ev -> muons.at(imu).pt     ;
      float glbptdiff = 0.;
      float trkptdiff = 0.;
      float qovertrkptdiff = 0.;
      float qoverglbptdiff = 0.;
      
      if ( ihltmu > -1){
	trkptdiff      = (mupttrk - ev->hltmuons.at(ihltmu).trkpt)/mupttrk;
	glbptdiff      = (muptglb - ev->hltmuons.at(ihltmu).pt)   /mupttrk;
	qovertrkptdiff = (1./mupttrk - 1./ev->hltmuons.at(ihltmu).trkpt)/(1./mupttrk);
	qoverglbptdiff = (1./muptglb - 1./ev->hltmuons.at(ihltmu).pt)   /(1./muptglb);
	
	PtHLT        ->Fill(trkptdiff     );
	QoverPtHLT   ->Fill(qovertrkptdiff);
	PtHLTGlb     ->Fill(glbptdiff     );
	QoverPtHLTGlb->Fill(qoverglbptdiff);
	
	if (fabs(mueta) <= 0.9 ){
	  PtHLT_barrel        ->Fill(trkptdiff     );
	  QoverPtHLT_barrel   ->Fill(qovertrkptdiff);
	  PtHLTGlb_barrel     ->Fill(glbptdiff     );
	  QoverPtHLTGlb_barrel->Fill(qoverglbptdiff);
	}
	if (fabs(mueta)>0.9 && fabs(mueta)<1.6){
	  PtHLT_int        ->Fill(trkptdiff     );
	  QoverPtHLT_int   ->Fill(qovertrkptdiff);
	  PtHLTGlb_int     ->Fill(glbptdiff     );
	  QoverPtHLTGlb_int->Fill(qoverglbptdiff);
	}
	if ( fabs(mueta)>=1.6){
	  PtHLT_endcap        ->Fill(trkptdiff     );
	  QoverPtHLT_endcap   ->Fill(qovertrkptdiff);
	  PtHLTGlb_endcap     ->Fill(glbptdiff     );
	  QoverPtHLTGlb_endcap->Fill(qoverglbptdiff);  
	}
      }	
      
      // Fill L3 muons from OI... 
      if ( itkmu > -1){
	trkptdiff      = (mupttrk - ev->tkmuons.at(itkmu).trkpt)/mupttrk;
	glbptdiff      = (muptglb - ev->tkmuons.at(itkmu).pt)   /mupttrk;
	qovertrkptdiff = (1./mupttrk - 1./ev->tkmuons.at(itkmu).trkpt)/(1./mupttrk);
	qoverglbptdiff = (1./muptglb - 1./ev->tkmuons.at(itkmu).pt)   /(1./muptglb);
	
	PtL3OI        ->Fill(trkptdiff     );
	QoverPtL3OI   ->Fill(qovertrkptdiff);
	PtL3OIGlb     ->Fill(glbptdiff     );
	QoverPtL3OIGlb->Fill(qoverglbptdiff);
	
	if (fabs(mueta) <= 0.9 ){
	  PtL3OI_barrel        ->Fill(trkptdiff     );
	  QoverPtL3OI_barrel   ->Fill(qovertrkptdiff);
	  PtL3OIGlb_barrel     ->Fill(glbptdiff     );
	  QoverPtL3OIGlb_barrel->Fill(qoverglbptdiff);
	}
	if (fabs(mueta)>0.9 && fabs(mueta)<1.6){
	  PtL3OI_int        ->Fill(trkptdiff     );
	  QoverPtL3OI_int   ->Fill(qovertrkptdiff);
	  PtL3OIGlb_int     ->Fill(glbptdiff     );
	  QoverPtL3OIGlb_int->Fill(qoverglbptdiff);
	}
	if ( fabs(mueta)>=1.6){
	  PtL3OI_endcap        ->Fill(trkptdiff     );
	  QoverPtL3OI_endcap   ->Fill(qovertrkptdiff);
	  PtL3OIGlb_endcap     ->Fill(glbptdiff     );
	  QoverPtL3OIGlb_endcap->Fill(qoverglbptdiff);  
	}
      }	

      // FILL HLT tracks for OI.
      if (ihltOI > -1){
	eta_muons_pass -> Fill (mueta);
	eta_hltOI_pass -> Fill (ev -> hltTrackOI.at(ihltOI).eta);
	eta_resta      -> Fill (mueta - ev -> hltTrackOI.at(ihltOI).eta);

	trkptdiff      = (mupttrk - ev->hltTrackOI.at(ihltOI).pt)/mupttrk;
	qovertrkptdiff = (1./mupttrk - 1./ev->hltTrackOI.at(ihltOI).pt)/(1./mupttrk);
	
	PtOI        ->Fill(trkptdiff     );
	QoverPtOI   ->Fill(qovertrkptdiff);
	
	if (fabs(mueta) <= 0.9 ){
	  eta_muons_barrel -> Fill (mueta);
	  eta_hltOI_pass_barrel -> Fill (ev -> hltTrackOI.at(ihltOI).eta);

	  PtOI_barrel        ->Fill(trkptdiff     );
	  QoverPtOI_barrel   ->Fill(qovertrkptdiff);
	}
	if (fabs(mueta)>0.9 && fabs(mueta)<1.6){
	  eta_muons_int -> Fill (mueta);
	  eta_hltOI_pass_int -> Fill (ev -> hltTrackOI.at(ihltOI).eta);

	  PtOI_int        ->Fill(trkptdiff     );
	  QoverPtOI_int   ->Fill(qovertrkptdiff);
	}
	if (fabs(mueta)>=1.6){
	  eta_muons_endcap -> Fill (mueta);
	  eta_hltOI_pass_endcap -> Fill (ev -> hltTrackOI.at(ihltOI).eta);

	  PtOI_endcap        ->Fill(trkptdiff     );
	  QoverPtOI_endcap   ->Fill(qovertrkptdiff);
	}
      }	
    }
  }

  // store histograms in a file.
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

  QoverPtHLT        -> Write();
  QoverPtHLT_barrel -> Write();
  QoverPtHLT_int    -> Write();
  QoverPtHLT_endcap -> Write();

  PtHLTGlb        -> Write();
  PtHLTGlb_barrel -> Write();
  PtHLTGlb_int    -> Write();
  PtHLTGlb_endcap -> Write();

  QoverPtHLTGlb        -> Write();
  QoverPtHLTGlb_barrel -> Write();
  QoverPtHLTGlb_int    -> Write();
  QoverPtHLTGlb_endcap -> Write();

  PtL3OI        -> Write();
  PtL3OI_barrel -> Write();
  PtL3OI_int    -> Write();
  PtL3OI_endcap -> Write();

  QoverPtL3OI        -> Write();
  QoverPtL3OI_barrel -> Write();
  QoverPtL3OI_int    -> Write();
  QoverPtL3OI_endcap -> Write();

  PtL3OIGlb        -> Write();
  PtL3OIGlb_barrel -> Write();
  PtL3OIGlb_int    -> Write();
  PtL3OIGlb_endcap -> Write();

  QoverPtL3OIGlb        -> Write();
  QoverPtL3OIGlb_barrel -> Write();
  QoverPtL3OIGlb_int    -> Write();
  QoverPtL3OIGlb_endcap -> Write();

  PtOI        -> Write();
  PtOI_barrel -> Write();
  PtOI_int    -> Write();
  PtOI_endcap -> Write();

  QoverPtOI        -> Write();
  QoverPtOI_barrel -> Write();
  QoverPtOI_int    -> Write();
  QoverPtOI_endcap -> Write();

  outfile    -> Close();  
  
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

void printProgBar( int percent ){
  std::string bar;  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}
