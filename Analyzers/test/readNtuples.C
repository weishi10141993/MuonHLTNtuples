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
#include "MuonHLT/Analyzers/src/MuonTree.h"
#include "TLorentzVector.h"

double muonmass = 0.10565837;

bool selectTagMuon  (MuonCand, TH1F* );
bool selectProbeMuon(MuonCand, MuonCand, TH1F* );
bool matchMuon      (MuonCand, std::vector<HLTObjCand>, std::string);
HLTMuonCand  matchL3        (MuonCand, std::vector<HLTMuonCand>);
L1MuonCand   matchL1        (MuonCand, std::vector<L1MuonCand>);

std::string L1filter      = "hltL1fL1sMu18L1Filtered0::TEST";
std::string L2filter      = "hltL2fL1sMu18L1f0L2Filtered10Q::TEST";
std::string L3filter      = "hltL3fL1sMu18L1f0L2f10QL3Filtered20Q::TEST";
std::string isofilter     = "hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09::TEST";
std::string isofilterTag  = "hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09::HLT";

double pt_bins[12]  = { 20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double iso_bins[12] = { 0  , 0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1   };

double offlineIsoCut = 0.15;

// ******************************************
//                                          *
//                                          *
int theRunNumber           = 274157;           
bool L2isfixed             = false;
std::string hltname        = "HLT_IsoMu20_v3";
std::string thepassfilter  = L3filter;
std::string theprobefilter = L1filter;
std::string effmeasured    = "L3";
//                                          *
//                                          *
// ******************************************

int miss_files[]       = {11, 13, 155, 156, 162, 171, 179, 193, 194, 196, 444, 458, 494, 198, 199, 243, 267, 301, 311, 383, 413, 420, 449, 472};
int miss_files_L2fix[] = {21, 32, 33, 35, 48, 51, 44, 67, 68, 82, 86, 99, 143, 147, 243, 241, 222, 201, 198};
int miss_files_273725[] = {134, 120, 191, 12, 14};

std::vector<int> missingfiles;
int firstfile, lastfile;
std::string folder, l2string;

void readNtuples(){

//   if (!L2isfixed) {
//     if (theRunNumber == 273725)  missingfiles.insert(missingfiles.end(), miss_files_273725, miss_files_273725+(sizeof(miss_files_273725)/sizeof(miss_files_273725[0])));
//     else                         missingfiles.insert(missingfiles.end(), miss_files, miss_files+(sizeof(miss_files)/sizeof(miss_files[0])));
//   }
//   else            missingfiles.insert(missingfiles.end(), miss_files_L2fix, miss_files_L2fix+(sizeof(miss_files_L2fix)/sizeof(miss_files_L2fix[0])));
    
  if (theRunNumber == 274157 && !L2isfixed){
    firstfile = 1;
    lastfile  = 2;
  } 
//   else if (theRunNumber == 273555){
//     firstfile = 136;
//     lastfile  = 153;
//   }
//   else if (theRunNumber == 273725){
//     firstfile = 1;
//     lastfile  = 50;
// //     lastfile  = 214;
//   }
//   else if (theRunNumber == 273555 && L2isfixed){
//     firstfile = 68;
//     lastfile  = 77;
//   } 
//   
//   if (!L2isfixed)  { 
//     if (theRunNumber == 273725) {folder = "crab_promptReco_v2_json_273725/160524_145139/"; l2string = "";}
//     else                        {folder = "crab_promptReco_v2_json_upTo273592/160523_090706/"; l2string = "";}
//   }
//   else             { folder = "crab_promptReco_v2_json_upTo273592_L2fix/160523_090802/"; l2string = "_L2fix";}
  
  folder = "../../../HLTrigger/Configuration/test/cfgs_2016_5_31_15_25/";
  bool doingL1 = thepassfilter.find("L1fL1") != std::string::npos;


  TFile* outfile = TFile::Open(Form("%sefficiency_%d%s.root", effmeasured.c_str(),theRunNumber,l2string.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;
  
  TH1F* dimuon_mass             = new TH1F("dimuon_mass"            ,"dimuon_mass"      , 1500,  0,  150 );
  TH1F* tagiso                  = new TH1F("tagiso"                 ,"tagiso"           ,  100,  0,  1   );
  TH1F* tagMuonPt               = new TH1F("tagMuonPt"              ,"tagMuonPt"        ,  150,  0,  150 );
  TH1F* nvtx_event              = new TH1F("nvtx_event"             ,"nvtx_event"       ,   60,  0,   60 );
 
  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"           ,   11,  pt_bins );
  TEfficiency* muonPt_barrel    = new TEfficiency("muonPt_barrel"   ,"muonPt_barrel"    ,   11,  pt_bins );
  TEfficiency* muonPt_endcap    = new TEfficiency("muonPt_endcap"   ,"muonPt_endcap"    ,   11,  pt_bins );
  TEfficiency* muonPt_eta0      = new TEfficiency("muonPt_eta0"     ,"muonPt_eta0"      ,   11,  pt_bins );
  TEfficiency* muonPt_eta1      = new TEfficiency("muonPt_eta1"     ,"muonPt_eta1"      ,   11,  pt_bins );
  TEfficiency* muonPt_eta2      = new TEfficiency("muonPt_eta2"     ,"muonPt_eta2"      ,   11,  pt_bins );
  
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"          ,   15, eta_bins );
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"          ,   40, -3.2, 3.2);
 
  TEfficiency* nvtx             = new TEfficiency("nvtx"            ,"nvtx"             ,   60,    0,  60);
  TEfficiency* nvtx_barrel      = new TEfficiency("nvtx_barrel"     ,"nvtx_barrel"      ,   60,    0,  60);
  TEfficiency* nvtx_endcap      = new TEfficiency("nvtx_endcap"     ,"nvtx_endcap"      ,   60,    0,  60);
  TEfficiency* nvtx_eta0        = new TEfficiency("nvtx_eta0"       ,"nvtx_eta0"        ,   60,    0,  60);
  TEfficiency* nvtx_eta1        = new TEfficiency("nvtx_eta1"       ,"nvtx_eta1"        ,   60,    0,  60);
  TEfficiency* nvtx_eta2        = new TEfficiency("nvtx_eta2"       ,"nvtx_eta2"        ,   60,    0,  60);
   
  TEfficiency* muonIso          = new TEfficiency("muonIso"         ,"muonIso"          ,  100,    0,   1);
  TEfficiency* muonIso_barrel   = new TEfficiency("muonIso_barrel"  ,"muonIso_barrel"   ,  100,    0,   1);
  TEfficiency* muonIso_endcap   = new TEfficiency("muonIso_endcap"  ,"muonIso_endcap"   ,  100,    0,   1);

  TH2F* muonPt_2d               = new TH2F("muonPt_2d"              ,"muonPt_2d"        ,   100,   0, 100, 100,    0, 100 );
  TH2F* muonEta_2d              = new TH2F("muonEta_2d"             ,"muonEta_2d"       ,   500, -2.5, 2.5,  500, -2.5, 2.5 );
//   TH2F* muonEta_2d              = new TH2F("muonEta_2d"             ,"muonEta_2d"       ,   50, -2.5, 2.5,  50, -2.5, 2.5 );
  TH2F* muonPhi_2d              = new TH2F("muonPhi_2d"             ,"muonPhi_2d"       ,   320,-3.2, 3.2, 320, -3.2, 3.2 );
   
  double offlineiso04 = 100;
  
  for (int ifile= firstfile; ifile < lastfile; ifile++){ 
    if ( std::find(missingfiles.begin(), missingfiles.end(), ifile) != missingfiles.end() ) continue;
    TFile* inputfile = TFile::Open(Form("%s/muonNtuples_274157_isolation.root",folder.c_str()),"READ");
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
  
    bool flagfile = false;
    for (Int_t eventNo=0; eventNo < nentries; eventNo++)
    {
      Int_t IgetEvent   = tree   -> GetEvent(eventNo);
      
      if (! (ev-> runNumber == theRunNumber)) break;
//       if (!flagfile){
//         std::cout << "input file with run 273555: " << inputfile -> GetName() << std::endl;
//         flagfile == true;
//       }
      unsigned int nmuons = ev->muons.size();
      if (nmuons < 2) continue;
      
      unsigned int nhltmuons = ev->hltmuons.size();
  //     if (nhltmuons > 0) std::cout << "Number of hlt muons = " << nhltmuons << std::endl;
      
      if (!ev-> hltTag.find(hltname)) continue;
      nvtx_event      -> Fill( ev -> nVtx   );
      
  
      for (int imu = 0; imu < nmuons; imu++){
        
        // select the tag muon        
        if (! selectTagMuon(ev -> muons.at(imu), tagiso))                            continue;
        if (! matchMuon(ev -> muons.at(imu), ev -> hltTag.objects, isofilterTag))    continue;
        tagMuonPt -> Fill(ev -> muons.at(imu).pt);
        
        for (int jmu = 0; jmu < nmuons; jmu++){
  
          bool pass   = false;
   
          // select the probe muon
          if (! selectProbeMuon(ev -> muons.at(jmu), ev -> muons.at(imu), dimuon_mass)) continue;
  
          offlineiso04 = ev -> muons.at(jmu).chargedDep_dR04 + std::max(0.,
                         ev -> muons.at(jmu).photonDep_dR04 + ev -> muons.at(jmu).neutralDep_dR04 - 0.5*ev -> muons.at(jmu).puPt_dR04);
          offlineiso04 = offlineiso04 / ev -> muons.at(jmu).pt;
          
          if (!(offlineiso04 < offlineIsoCut)) continue;


          if (!doingL1 && !(matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, theprobefilter))) continue;
  
          // match probe muon to the interesting filter and fill numerator histograms
          if (matchMuon(ev -> muons.at(jmu), ev -> hlt.objects, thepassfilter))       pass = true; 
  
          if (pass){
//             L1MuonCand theL3 = matchL1(ev -> muons.at(jmu), ev -> L1muons);
            HLTMuonCand theL3 = matchL3(ev -> muons.at(jmu), ev -> L2muons);
//             HLTMuonCand theL3 = matchL3(ev -> muons.at(jmu), ev -> hltmuons);
            muonPt_2d  -> Fill( ev -> muons.at(jmu).pt,  theL3.pt );
            muonEta_2d -> Fill( ev -> muons.at(jmu).eta, theL3.eta);
            muonPhi_2d -> Fill( ev -> muons.at(jmu).phi, theL3.phi);
          }
          
          muonPt       -> Fill( pass, ev -> muons.at(jmu).pt );
          muonEta      -> Fill( pass, ev -> muons.at(jmu).eta);
          muonPhi      -> Fill( pass, ev -> muons.at(jmu).phi);
          muonIso      -> Fill( pass, offlineiso04           );
          nvtx         -> Fill( pass, ev -> nVtx             );
          
          if (fabs(ev -> muons.at(jmu).eta) < 1.479){
            muonPt_barrel     -> Fill( pass, ev -> muons.at(jmu).pt );
            nvtx_barrel       -> Fill( pass, ev -> nVtx );
            muonIso_barrel    -> Fill( pass, offlineiso04 );
          }
          else{
            muonPt_endcap     -> Fill( pass, ev -> muons.at(jmu).pt );
            nvtx_endcap       -> Fill( pass, ev -> nVtx             );
            muonIso_endcap    -> Fill( pass, offlineiso04           );
          }
  
//           if (fabs(ev -> muons.at(jmu).eta) < 0.9){
//             muonPt_eta0     -> Fill( pass, ev -> muons.at(jmu).pt );
//             nvtx_eta0       -> Fill( pass, ev -> nVtx );
//           }
//           else if (fabs(ev -> muons.at(jmu).eta) > 0.9 && fabs(ev -> muons.at(jmu).eta) < 1.2){
//             muonPt_eta1     -> Fill( pass, ev -> muons.at(jmu).pt );
//             nvtx_eta1       -> Fill( pass, ev -> nVtx );
//           }
//           else{
//             muonPt_eta2     -> Fill( pass, ev -> muons.at(jmu).pt );
//             nvtx_eta2       -> Fill( pass, ev -> nVtx );
//           }
        }
      }
    }  
  } //end loop files
  
  outfile           -> cd();
  muonPt            -> Write();
  muonPt_barrel     -> Write();
  muonPt_endcap     -> Write();
  muonPt_eta0       -> Write();
  muonPt_eta1       -> Write();
  muonPt_eta2       -> Write();
  tagMuonPt         -> Write();
  
  muonEta           -> Write();
  muonPhi           -> Write();
  
  nvtx_event        -> Write();
  nvtx              -> Write();
  nvtx_barrel       -> Write();
  nvtx_endcap       -> Write();
  nvtx_eta0         -> Write();
  nvtx_eta1         -> Write();
  nvtx_eta2         -> Write();
  
  muonIso           -> Write();
  muonIso_barrel    -> Write();
  muonIso_endcap    -> Write();

  dimuon_mass       -> Write();
  tagiso            -> Write();
  
  muonPt_2d         -> Write();
  muonEta_2d        -> Write();
  muonPhi_2d        -> Write();

  
  outfile           -> Close();  
  
  return;
}



bool matchMuon(MuonCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

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


bool selectTagMuon(MuonCand mu, TH1F* tagiso){
  
  if (!( mu.pt         > 22  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  //add isolation cut
  float offlineiso04 = mu.chargedDep_dR04 + std::max(0.,
                       mu.photonDep_dR04 + mu.neutralDep_dR04 - 0.5*mu.puPt_dR04);
  offlineiso04       = offlineiso04 / mu.pt;
  tagiso -> Fill(offlineiso04);
  if (!(offlineiso04   < offlineIsoCut)) return false;
  
  return true;
}


bool selectProbeMuon(MuonCand mu, MuonCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt == tagMu.pt  && 
      mu.pt == tagMu.eta &&
      mu.pt == tagMu.phi ) 
    return false;
  
  if (!( mu.pt         > 22  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (!( mu.isTight    == 1  )) return false; 
  
  if (mu.charge * tagMu.charge > 0) return false;
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();
  dimuon_mass -> Fill(mumumass);
  if (! (mumumass > 86. && mumumass < 96. )) return false;
  
  return true;
}



HLTMuonCand matchL3(MuonCand mu, std::vector<HLTMuonCand> L3cands){

  bool match = false;
  int nL3 = L3cands.size();

  float minDR = 0.1;
  float theDR = 100;
  HLTMuonCand theL3;
  theL3.pt        = -1000;
  theL3.eta       = -1000;
  theL3.phi       = -1000;
  theL3.trkpt     = -1000;
  theL3.ecalDep   = -1000;
  theL3.hcalDep   = -1000;
  theL3.trkDep    = -1000;
  theL3.ecalDep05 = -1000;
  theL3.hcalDep05 = -1000;
  theL3.ecalDep1  = -1000;
  theL3.hcalDep1  = -1000;
  
  for ( std::vector<HLTMuonCand>::const_iterator it = L3cands.begin(); it != L3cands.end(); ++it ) {
      
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
      theL3 = *it;
    }
  }
  
  return theL3;
}


L1MuonCand matchL1(MuonCand mu, std::vector<L1MuonCand> L1cands){

  bool match = false;
  int nL1 = L1cands.size();

  float minDR = 0.3;
  float theDR = 100;
  L1MuonCand theL1;
  theL1.pt        = -1000;
  theL1.eta       = -1000;
  theL1.phi       = -1000;
  
  for ( std::vector<L1MuonCand>::const_iterator it = L1cands.begin(); it != L1cands.end(); ++it ) {
      
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
    if (theDR < minDR){
      minDR = theDR;
      match = true;
      theL1 = *it;
    }
  }
  
  return theL1;
}


