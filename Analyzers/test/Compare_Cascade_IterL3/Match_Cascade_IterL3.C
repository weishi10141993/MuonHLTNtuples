
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

vector<size_t> matchHltTrackMu ( std::vector<HltTrackCand>, std::vector<HltTrackCand> );
bool matchMuon (HltTrackCand , std::vector<HLTObjCand>, std::string);

std::string thepassfilter = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 

void Match_Cascade_IterL3 (TString file_cascade = "/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/old_reco/muonNtuples.root", TString file_iterl3 = "/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/newIterL3/muonNtuples.root") {
  
  TFile* outfile = TFile::Open("FilterL3_Cascade_IterL3_Match_NOTHG.root", "RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  double pt_bins[17]  = { 5, 7, 9, 12, 16,  20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
  double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};

  //Create histograms   
  TH1F*eta_position = new TH1F("h_eta_position", "eta_position", 3, -0.5, 2.5);
  TH1F* eta_cas   = new TH1F("h_eta_cas"   , "eta_cas"   ,   14,  eta_bins); 
  TH1F* Pt_cas     = new TH1F("h_Pt_cas"     , "Pt_cas"     ,   16,  pt_bins ); 
  TH1F* phi_cas    = new TH1F("h_phi_cas"    , "phi_cas"    ,   20,  -3.2, 3.2);
  TH1F* dxy_cas    = new TH1F("h_dxy_cas"    , "dxy_cas"    ,  100,  -0.3, 0.3);
  TH1F* dz_cas     = new TH1F("h_dz_cas"     , "dz_cas"     ,   20,   -3.,  3.);
  TH1F* chi2_cas   = new TH1F("h_chi2_cas"   , "chi2_cas"   ,   60,    0.,  7.);
  TH1F* PixHit_cas = new TH1F("h_PixHit_cas" , "PixHit_cas" ,   20,  -0.5,19.5);
  TH1F* LayHit_cas = new TH1F("h_LayHit_cas" , "LayHit_cas" ,   16,   2.5,18.5);
  TH1F* PixLay_cas = new TH1F("h_PixLay_cas" , "PixLay_cas" ,    7,  -0.5, 6.5);
  TH1F* FracValTrackHit_cas = new TH1F("h_FracValTrackHit_cas" , "FracValTrackHit_cas" ,    100,  0.5, 1.5);

  TH1F* eta_delta   = new TH1F("h_eta_delta"   , "eta_delta"   ,   14,  eta_bins); 
  TH1F* Pt_delta     = new TH1F("h_Pt_delta"     , "Pt_delta"     ,   100, -10.0, 10.0 ); 
  TH1F* phi_delta    = new TH1F("h_phi_delta"    , "phi_delta"    ,   20,  -3.2, 3.2);
  TH1F* dxy_delta    = new TH1F("h_dxy_delta"    , "dxy_delta"    ,  100,  -0.03, 0.03);
  TH1F* dz_delta     = new TH1F("h_dz_delta"     , "dz_delta"     ,   50,   -0.2,  0.2);
  TH1F* chi2_delta   = new TH1F("h_chi2_delta"   , "chi2_delta"   ,   80,    -1.5,  1.5);
  TH1F* PixHit_delta = new TH1F("h_PixHit_delta" , "PixHit_delta" ,   5,  -2.5,2.5);
  TH1F* LayHit_delta = new TH1F("h_LayHit_delta" , "LayHit_delta" ,   11,  - 5.5,5.5);
  TH1F* PixLay_delta = new TH1F("h_PixLay_delta" , "PixLay_delta" ,    5,  -2.5, 2.5);
  TH1F* FracValTrackHit_delta = new TH1F("h_FracValTrackHit_delta" , "FracValTrackHit_delta" ,    100,  -0.5, 0.5);
  
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
 
    for (Int_t event_cas=0; event_cas < nentries_cas ;  event_cas++) {
      Int_t getEvent_cas = tree_cas -> GetEvent(event_cas);

      if ((ev_cas -> eventNumber) != (ev_iter -> eventNumber)) continue;
      if ((ev_cas -> runNumber)   != (ev_iter -> runNumber))   continue;
      if ((ev_cas -> luminosityBlockNumber) != (ev_iter -> luminosityBlockNumber)) continue;

      unsigned int nhltOI_cas   = ev_cas -> hltTrackOI.size(); 
      vector< size_t> matchTrackOI   = matchHltTrackMu(ev_cas -> hltTrackOI, ev_iter -> hltTrackOI );

      for (size_t ihlt_cas = 0; ihlt_cas < nhltOI_cas; ihlt_cas++){
	if (!matchMuon(ev_cas -> hltTrackOI.at(ihlt_cas), ev_cas -> hlt.objects, thepassfilter))  continue;
	int matching = matchTrackOI [ihlt_cas];
	  	  
	if ( matching < 0) {
   	  eta_cas   -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).eta );
	  phi_cas    -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).phi );
	  Pt_cas     -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).pt );
	  chi2_cas   -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).chi2 );
	  dxy_cas    -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).dxy );
	  dz_cas     -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).dz );
	  PixHit_cas -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).pixelHits );
	  LayHit_cas -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).layerHits );
	  PixLay_cas -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).pixelLayers);
	  FracValTrackHit_cas -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).fracValidTrackhit );
	}

	if (nhltOI_cas != nhltOI_iter && matching > -1){
	  if (!matchMuon(ev_iter -> hltTrackOI.at(matching), ev_iter -> hlt.objects, thepassfilter)) continue;
   	  eta_delta   -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).eta         - ev_iter -> hltTrackOI.at(matching).eta);
	  phi_delta    -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).phi         - ev_iter -> hltTrackOI.at(matching).phi );
	  Pt_delta     -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).pt          - ev_iter -> hltTrackOI.at(matching).pt );
	  chi2_delta   -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).chi2        - ev_iter -> hltTrackOI.at(matching).chi2 );
	  dxy_delta    -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).dxy         - ev_iter -> hltTrackOI.at(matching).dxy);
	  dz_delta     -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).dz          - ev_iter -> hltTrackOI.at(matching).dz);
	  PixHit_delta -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).pixelHits   - ev_iter -> hltTrackOI.at(matching).pixelHits);
	  LayHit_delta -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).layerHits   - ev_iter -> hltTrackOI.at(matching).layerHits);
	  PixLay_delta -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).pixelLayers - ev_iter -> hltTrackOI.at(matching).pixelLayers);
	  FracValTrackHit_delta -> Fill ( ev_cas -> hltTrackOI.at(ihlt_cas).fracValidTrackhit - ev_iter -> hltTrackOI.at(matching).fracValidTrackhit);

	}
      }
    }
  }
  outfile -> cd();

  eta_position -> Write();

  eta_cas   -> Write();
  phi_cas    -> Write();
  Pt_cas     -> Write();
  chi2_cas   -> Write();
  dxy_cas    -> Write();
  dz_cas     -> Write();
  PixHit_cas -> Write();
  LayHit_cas -> Write();
  PixLay_cas -> Write();
  FracValTrackHit_cas -> Write ();


  eta_delta   -> Write();
  phi_delta    -> Write();
  Pt_delta     -> Write();
  chi2_delta   -> Write();
  dxy_delta    -> Write();
  dz_delta     -> Write();
  PixHit_delta -> Write();
  LayHit_delta -> Write();
  PixLay_delta -> Write();
  FracValTrackHit_delta -> Write ();
 
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




vector<size_t> matchHltTrackMu (  std::vector<HltTrackCand> hltTrack1, std::vector<HltTrackCand> hltTrack2 ){

  unsigned int nmuons   = hltTrack1.size();
  unsigned int nhlt     = hltTrack2.size(); 
 
  float deltaRMat[nmuons][nhlt];
  double theDeltaR;
  int mu_counter = -2;
  int hlt_counter = -2;

  for (int imu = 0; imu < nmuons; imu++){ 
    for ( int ihlt=0; ihlt<hltTrack2.size(); ihlt++ ) {

      deltaRMat[imu][ihlt]=deltaR( hltTrack1.at(imu).eta, hltTrack1.at(imu).phi, hltTrack2.at(ihlt).eta, hltTrack2.at(ihlt).phi);
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
      for ( int ihlt=0; ihlt < hltTrack2.size(); ihlt++ ) {
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

