//Counts numbers of HLT muon from Cascade that:   0.No match with offline muon   1.Match HLTmuon  2.Match TkMuon  3. Tk or HLTmuon  4.Passes HLT filter 5.Passes Tk filter    6.Or in the filter.

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

double muonmass = 0.10565837;

bool selectTagMuon  (MuonCand);
bool selectProbeMuon(MuonCand, MuonCand);
bool selectMuon     (MuonCand);
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
bool  matchMuonWithL3 (MuonCand, std::vector<HLTMuonCand>);
void printProgBar(int);

std::string L1filter      =  "hltL1fL1sMu22or25L1Filtered0::TEST"; 
std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
std::string L3filter      =  "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 
std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07::HLT";

double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double offlineIsoCut = 0.15;

// ******************************************
//                                          *
//                                          *
std::string hltname        = "HLT_IsoMu27_v10";
std::string thepassfilter  = L3filter;
std::string thepassfilter2 = "hltL3fL1sMu22Or25f0TkFiltered27Q::TEST";
std::string theprobefilter = L1filter; 
float offlinePtCut         = 30.;
//                                          *
//                                          *
// ******************************************

void testCascade(TString inputfilename="/afs/cern.ch/user/s/sferrere/private/CMSSW_9_2_7/src/workspace/old_reco/results/result_old_reco.root", std::string effmeasured="test"){

  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos; 

  TFile* outfile = TFile::Open(Form("%s_Cascade.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  TH1F* counter    = new TH1F("h_counter" ,"counter" ,   7,  -0.5,   6.5  );
  TH1F* eta        = new TH1F("h_eta"     ,"eta"     ,  15,  eta_bins);
 
  double offlineiso04 = 100;
  
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

  for (Int_t eventNo=0; eventNo < nentries; eventNo++) 
    {
      Int_t IgetEvent   = tree   -> GetEvent(eventNo); 
      printProgBar((int)(eventNo*100./nentries));
    
      unsigned int nmuons = ev->muons.size(); 
      if (nmuons < 2) continue; 
      unsigned int nhltmuons = ev->hltmuons.size(); 
    
      for (int imu = 0; imu < nmuons; imu++){ 
        
	if (! selectTagMuon(ev -> muons.at(imu))) continue; 

	if (! matchMuon(ev -> muons.at(imu), ev -> hltTag.objects, isofilterTag)) continue;

	for (int jmu = imu+1; jmu < nmuons; jmu++){
	  bool pass   = false;

	  // select the probe muon
	  if (!selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu))) continue;

	  if (!doingL1 && !(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, theprobefilter))) continue;
	  
	  //FILL COUNTER HISTOGRAM
	  int var=0;
	  //Prefilter part
	  if ( matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons)){ 
	    var = 1; 
	    counter -> Fill ( var );
	  }

	  if ( matchMuonWithL3(ev->muons.at(jmu),ev->tkmuons)) {
	    var = 2; 
	    counter -> Fill ( var );
	  }

	  if ( matchMuonWithL3(ev->muons.at(jmu),ev->hltmuons) || matchMuonWithL3(ev->muons.at(jmu),ev->tkmuons)) {
	    var = 3; 
	    counter -> Fill ( var );
	  }


	  //POST FILTER
	  if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, thepassfilter)){ 
	    var = 4; 
	    counter -> Fill ( var );
	  }

	  if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, thepassfilter2)){ 
	    var = 5; 
	    counter -> Fill ( var );
	  }

	  if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, thepassfilter) || matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, thepassfilter2)){ 
	    var = 6; 
	    counter -> Fill ( var );
	  }

	  //NO MATCH
	  if ( var == 0 ) counter -> Fill ( var );

	} 
      }
    }  

  outfile -> cd();

  eta     -> Write();
  counter -> Write();

  outfile -> Close();  
  
  return;
}


bool matchMuon(MuonCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1;
  if (tagFilterName.find("L1fL1") != std::string::npos) minDR = 1.0;
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



bool selectTagMuon(MuonCand mu){
  
  if (!( mu.pt         > 24  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0., mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt; //calculate isolation value

  if (offlineiso04   > offlineIsoCut) return false; //if the isolation requirement is not satisfied is not a good tag muon

  return true;
}




bool selectMuon(MuonCand mu){  
  if (!( mu.pt         > offlinePtCut  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false;
  if (!( mu.isLoose    == 1  )) return false; 
  return true;
}




bool selectProbeMuon(MuonCand mu, MuonCand tagMu){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt          > 0  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  if (mu.charge * tagMu.charge > 0) return false;

  float offlineiso04 = mu.chargedDep_dR04 + std::max(0., mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  if (offlineiso04   > offlineIsoCut) return false; 
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();

  if (! (mumumass > 86. && mumumass < 96. )) return false;
  
  return true;
}



bool matchMuonWithL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

  bool match = false;
  float minDR = 0.3;
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
