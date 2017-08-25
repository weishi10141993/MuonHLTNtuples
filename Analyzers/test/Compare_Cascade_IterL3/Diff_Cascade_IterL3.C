
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

vector<size_t> matchHltTrackMu ( std::vector<MuonCand>, std::vector<HltTrackCand> );
bool matchMuon (HltTrackCand , std::vector<HLTObjCand>, std::string);

//std::string thepassfilter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 

void Diff_Cascade_IterL3 (TString file_cascade = "/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/old_reco/muonNtuples.root", TString file_iterl3 = "/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/newIterL3/oldJSON/muonNtuples_HP.root") {
  
  TFile* outfile = TFile::Open("TrackDiff_Cascade_IterL3.root", "RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  double pt_bins[17]  = { 5, 7, 9, 12, 16,  20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
  double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};

  //Create histograms   
  TH1F* eta_pass   = new TH1F("h_eta_pass"   , "eta_pass"   ,   14,  eta_bins); 
  TH1F* Pt_pass     = new TH1F("h_Pt_pass"     , "Pt_pass"     ,   16,  pt_bins ); 
  TH1F* phi_pass    = new TH1F("h_phi_pass"    , "phi_pass"    ,   20,  -3.2, 3.2);
  TH1F* dxy_pass    = new TH1F("h_dxy_pass"    , "dxy_pass"    ,  100,  -0.3, 0.3);
  TH1F* dz_pass     = new TH1F("h_dz_pass"     , "dz_pass"     ,   20,   -3.,  3.);
  TH1F* chi2_pass   = new TH1F("h_chi2_pass"   , "chi2_pass"   ,   60,    0.,  7.);
  TH1F* PixHit_pass = new TH1F("h_PixHit_pass" , "PixHit_pass" ,   20,  -0.5,19.5);
  TH1F* LayHit_pass = new TH1F("h_LayHit_pass" , "LayHit_pass" ,   16,   2.5,18.5);
  TH1F* PixLay_pass = new TH1F("h_PixLay_pass" , "PixLay_pass" ,    7,  -0.5, 6.5);
  TH1F* FracValTrackHit_pass = new TH1F("h_FracValTrackHit_pass" , "FracValTrackHit_pass" ,    100,  0.5, 1.5);

  TH1F* eta_pass_barrel   = new TH1F("h_eta_pass_barrel"   , "eta_pass_barrel"   ,   14,  eta_bins); 
  TH1F* Pt_pass_barrel     = new TH1F("h_Pt_pass_barrel"     , "Pt_pass_barrel"     ,   16,  pt_bins ); 
  TH1F* phi_pass_barrel    = new TH1F("h_phi_pass_barrel"    , "phi_pass_barrel"    ,   20,  -3.2, 3.2);
  TH1F* dxy_pass_barrel    = new TH1F("h_dxy_pass_barrel"    , "dxy_pass_barrel"    ,  100,  -0.3, 0.3);
  TH1F* dz_pass_barrel     = new TH1F("h_dz_pass_barrel"     , "dz_pass_barrel"     ,   20,   -3.,  3.);
  TH1F* chi2_pass_barrel   = new TH1F("h_chi2_pass_barrel"   , "chi2_pass_barrel"   ,   60,    0.,  7.);
  TH1F* PixHit_pass_barrel = new TH1F("h_PixHit_pass_barrel" , "PixHit_pass_barrel" ,   20,  -0.5,19.5);
  TH1F* LayHit_pass_barrel = new TH1F("h_LayHit_pass_barrel" , "LayHit_pass_barrel" ,   16,   2.5,18.5);
  TH1F* PixLay_pass_barrel = new TH1F("h_PixLay_pass_barrel" , "PixLay_pass_barrel" ,    7,  -0.5, 6.5);
  TH1F* FracValTrackHit_pass_barrel = new TH1F("h_FracValTrackHit_pass_barrel" , "FracValTrackHit_pass_barrel" ,    100,  0.5, 1.5);

  TH1F* eta_pass_int   = new TH1F("h_eta_pass_int"   , "eta_pass_int"   ,   14,  eta_bins); 
  TH1F* Pt_pass_int     = new TH1F("h_Pt_pass_int"     , "Pt_pass_int"     ,   16,  pt_bins ); 
  TH1F* phi_pass_int    = new TH1F("h_phi_pass_int"    , "phi_pass_int"    ,   20,  -3.2, 3.2);
  TH1F* dxy_pass_int    = new TH1F("h_dxy_pass_int"    , "dxy_pass_int"    ,  100,  -0.3, 0.3);
  TH1F* dz_pass_int     = new TH1F("h_dz_pass_int"     , "dz_pass_int"     ,   20,   -3.,  3.);
  TH1F* chi2_pass_int   = new TH1F("h_chi2_pass_int"   , "chi2_pass_int"   ,   60,    0.,  7.);
  TH1F* PixHit_pass_int = new TH1F("h_PixHit_pass_int" , "PixHit_pass_int" ,   20,  -0.5,19.5);
  TH1F* LayHit_pass_int = new TH1F("h_LayHit_pass_int" , "LayHit_pass_int" ,   16,   2.5,18.5);
  TH1F* PixLay_pass_int = new TH1F("h_PixLay_pass_int" , "PixLay_pass_int" ,    7,  -0.5, 6.5);
  TH1F* FracValTrackHit_pass_int = new TH1F("h_FracValTrackHit_pass_int" , "FracValTrackHit_pass_int" ,    100,  0.5, 1.5);

  TH1F* eta_pass_endcap   = new TH1F("h_eta_pass_endcap"   , "eta_pass_endcap"   ,   14,  eta_bins); 
  TH1F* Pt_pass_endcap     = new TH1F("h_Pt_pass_endcap"     , "Pt_pass_endcap"     ,   16,  pt_bins ); 
  TH1F* phi_pass_endcap    = new TH1F("h_phi_pass_endcap"    , "phi_pass_endcap"    ,   20,  -3.2, 3.2);
  TH1F* dxy_pass_endcap    = new TH1F("h_dxy_pass_endcap"    , "dxy_pass_endcap"    ,  100,  -0.3, 0.3);
  TH1F* dz_pass_endcap     = new TH1F("h_dz_pass_endcap"     , "dz_pass_endcap"     ,   20,   -3.,  3.);
  TH1F* chi2_pass_endcap   = new TH1F("h_chi2_pass_endcap"   , "chi2_pass_endcap"   ,   60,    0.,  7.);
  TH1F* PixHit_pass_endcap = new TH1F("h_PixHit_pass_endcap" , "PixHit_pass_endcap" ,   20,  -0.5,19.5);
  TH1F* LayHit_pass_endcap = new TH1F("h_LayHit_pass_endcap" , "LayHit_pass_endcap" ,   16,   2.5,18.5);
  TH1F* PixLay_pass_endcap = new TH1F("h_PixLay_pass_endcap" , "PixLay_pass_endcap" ,    7,  -0.5, 6.5);
  TH1F* FracValTrackHit_pass_endcap = new TH1F("h_FracValTrackHit_pass_endcap" , "FracValTrackHit_pass_endcap" ,    100,  0.5, 1.5);

  
  TFile* input_cascade = TFile::Open(file_cascade, "READ"); 
  std::cout << "input file Cascade: " << input_cascade  -> GetName() << std::endl;
  TFile* input_iterl3 = TFile::Open(file_iterl3, "READ"); 
  std::cout << "input file Iterative L3: " << input_iterl3 -> GetName() << std::endl;

  TTree *tree_iter = (TTree*) input_iterl3  -> Get("muonNtuples/muonTree"); 
  TTree *tree_cas  = (TTree*) input_cascade -> Get("muonNtuples/muonTree"); 
    
  if (!tree_iter || !tree_cas) {
    std::cout << " *** tree not found *** " << std::endl;
    return;
  }
    
  MuonEvent* ev_iter      = new MuonEvent(); 
  TBranch*  evBranch_iter = tree_iter -> GetBranch("event");
  evBranch_iter -> SetAddress(&ev_iter); 

  MuonEvent* ev_cas      = new MuonEvent(); 
  TBranch*  evBranch_cas = tree_cas -> GetBranch("event");
  evBranch_cas -> SetAddress(&ev_cas);

  int nentries_cas   = tree_cas -> GetEntriesFast();
  std::cout << "Number of entries of Cascade+TkMu = " << nentries_cas << std::endl;
  int nentries_iter = tree_iter -> GetEntriesFast();
  std::cout << "Number of entries of iterative L3 = " << nentries_iter << std::endl;


  for (Int_t event_iter=0; event_iter < nentries_iter; event_iter++) {
    Int_t getEvent_iter = tree_iter -> GetEvent(event_iter);
    unsigned int nhltOI_iter   = ev_iter -> hltTrackOI.size();
    unsigned int nmuons_iter   = ev_iter -> muons.size();

 
    for (Int_t event_cas=0; event_cas < nentries_cas ;  event_cas++) {
      Int_t getEvent_cas = tree_cas -> GetEvent(event_cas);

      if ((ev_cas -> eventNumber) != (ev_iter -> eventNumber)) continue;
      if ((ev_cas -> runNumber)   != (ev_iter -> runNumber))   continue;
      if ((ev_cas -> luminosityBlockNumber) != (ev_iter -> luminosityBlockNumber)) continue;


      for (size_t imuon_iter = 0; imuon_iter < nmuons_iter; imuon_iter++){

	//if (!matchMuon(ev_cas -> hltTrackOI.at(ihlt_cas), ev_cas -> hlt.objects, thepassfilter))  continue;

	vector< size_t> matchTrackOI_iter = matchHltTrackMu(ev_iter -> muons, ev_iter -> hltTrackOI );
	int matching_iter = matchTrackOI_iter[imuon_iter];

	vector< size_t> matchTrackOI_cas   = matchHltTrackMu(ev_iter -> muons, ev_cas -> hltTrackOI );
	int matching_cas = matchTrackOI_cas[imuon_iter];

	if (matching_iter > -1) continue;

	if (matching_iter < 0 ){
	  
	  if (matching_cas > -1) {
	    eta_pass    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).eta );
	    phi_pass    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).phi );
	    Pt_pass     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pt );
	    chi2_pass   -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).chi2 );
	    dxy_pass    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dxy );
	    dz_pass     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dz );
	    PixHit_pass -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelHits );
	    LayHit_pass -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).layerHits );
	    PixLay_pass -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelLayers);
	    FracValTrackHit_pass -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).fracValidTrackhit );
	    
	    if ( fabs(ev_iter -> muons.at(imuon_iter).eta) <= 0.9 ) {
	      eta_pass_barrel    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).eta );
	      phi_pass_barrel    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).phi );
	      Pt_pass_barrel     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pt );
	      chi2_pass_barrel   -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).chi2 );
	      dxy_pass_barrel    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dxy );
	      dz_pass_barrel     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dz );
	      PixHit_pass_barrel -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelHits );
	      LayHit_pass_barrel -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).layerHits );
	      PixLay_pass_barrel -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelLayers);
	      FracValTrackHit_pass_barrel -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).fracValidTrackhit );
	    }
	    
	    if ( fabs(ev_iter -> muons.at(imuon_iter).eta) > 0.9 && fabs(ev_iter -> muons.at(imuon_iter).eta) < 1.6) {
	      eta_pass_int    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).eta );
	      phi_pass_int    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).phi );
	      Pt_pass_int     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pt );
	      chi2_pass_int   -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).chi2 );
	      dxy_pass_int    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dxy );
	      dz_pass_int     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dz );
	      PixHit_pass_int -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelHits );
	      LayHit_pass_int -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).layerHits );
	      PixLay_pass_int -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelLayers);
	      FracValTrackHit_pass_int -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).fracValidTrackhit );
	    }
	    
	    if (fabs(ev_iter -> muons.at(imuon_iter).eta) >= 1.6) {
	      eta_pass_endcap    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).eta );
	      phi_pass_endcap    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).phi );
	      Pt_pass_endcap     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pt );
	      chi2_pass_endcap   -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).chi2 );
	      dxy_pass_endcap    -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dxy );
	      dz_pass_endcap     -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).dz );
	      PixHit_pass_endcap -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelHits );
	      LayHit_pass_endcap -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).layerHits );
	      PixLay_pass_endcap -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).pixelLayers);
	      FracValTrackHit_pass_endcap -> Fill ( ev_cas -> hltTrackOI.at(matching_cas).fracValidTrackhit );
	    }
	  }
	}
      }
    }
  }
  outfile -> cd();

  eta_pass   -> Write();
  phi_pass    -> Write();
  Pt_pass     -> Write();
  chi2_pass   -> Write();
  dxy_pass    -> Write();
  dz_pass     -> Write();
  PixHit_pass -> Write();
  LayHit_pass -> Write();
  PixLay_pass -> Write();
  FracValTrackHit_pass -> Write ();

  eta_pass_barrel   -> Write();
  phi_pass_barrel    -> Write();
  Pt_pass_barrel     -> Write();
  chi2_pass_barrel   -> Write();
  dxy_pass_barrel    -> Write();
  dz_pass_barrel     -> Write();
  PixHit_pass_barrel -> Write();
  LayHit_pass_barrel -> Write();
  PixLay_pass_barrel -> Write();
  FracValTrackHit_pass_barrel -> Write ();

  eta_pass_int   -> Write();
  phi_pass_int    -> Write();
  Pt_pass_int     -> Write();
  chi2_pass_int   -> Write();
  dxy_pass_int    -> Write();
  dz_pass_int     -> Write();
  PixHit_pass_int -> Write();
  LayHit_pass_int -> Write();
  PixLay_pass_int -> Write();
  FracValTrackHit_pass_int -> Write ();

  eta_pass_endcap   -> Write();
  phi_pass_endcap    -> Write();
  Pt_pass_endcap     -> Write();
  chi2_pass_endcap   -> Write();
  dxy_pass_endcap    -> Write();
  dz_pass_endcap     -> Write();
  PixHit_pass_endcap -> Write();
  LayHit_pass_endcap -> Write();
  PixLay_pass_endcap -> Write();
  FracValTrackHit_pass_endcap -> Write ();
 
  outfile      -> Close();  
  
  return;

}



bool matchMuon(HltTrackCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1;
  if (tagFilterName.find("L1fL1") != std::string::npos) minDR = 0.3;
  float theDR = 100;

  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) {

    if ( it->filterTag.compare(tagFilterName) == 0) { 
      theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);

      if (theDR < minDR){
        minDR = theDR;
        match = true;
      }
    }
  }
  
  return match;
}




vector<size_t> matchHltTrackMu (  std::vector<MuonCand> mu, std::vector<HltTrackCand> hltTrack ){

  unsigned int nmuons   = mu.size();
  unsigned int nhlt     = hltTrack.size(); 
 
  float deltaRMat[nmuons][nhlt];
  double theDeltaR;
  int mu_counter = -2;
  int hlt_counter = -2;

  for (int imu = 0; imu < nmuons; imu++){ 
    for ( int ihlt=0; ihlt<hltTrack.size(); ihlt++ ) {

      deltaRMat[imu][ihlt]=deltaR( mu.at(imu).eta, mu.at(imu).phi, hltTrack.at(ihlt).eta, hltTrack.at(ihlt).phi);
      theDeltaR = deltaRMat[imu][ihlt];
    }
  }

  float maxDeltaR = 0.1;
  std::vector<size_t>  match(nmuons, -1);

  for (size_t k=0; k < nmuons; k++){
    double minDeltaR = maxDeltaR;
    size_t muon = -1;
    size_t hlt = -1;

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

