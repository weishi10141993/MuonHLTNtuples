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
bool debug = false;

enum Sig { 
  Prompt = 0,
  DiMuon,
  LowPt,
  DisplacedOld,
  DisplacedNew,
};

bool selectTagMuon      (GenParticleCand, TH1F* );
bool selectProbeMuon    (GenParticleCand, GenParticleCand, TH1F* );
bool selectMuon         (GenParticleCand);
bool matchMuon          (GenParticleCand, std::vector<HLTObjCand>, std::string);
bool matchMuonWithHLT   (GenParticleCand, std::vector<HLTMuonCand>);
bool matchMuonWithL1    (GenParticleCand, std::vector<L1MuonCand>);
bool matchHLTMuonWithGen(HLTMuonCand, std::vector<GenParticleCand>);
bool matchL1MuonWithGen (L1MuonCand,  std::vector<GenParticleCand>);
float getLeadingPtCut(int);
float getTrailingPtCut(int);
void printProgBar(int);

double pt_bins[17]  = { 5, 7, 9, 12, 16,  20 ,  24 ,  27 ,   30,   35,   40,   45,   50,  60, 70 ,  90, 150 };
double dz_bins[11]  = {-15, -8, -6, -4, -2, 0, 2, 4, 6, 8, 15};
double eta_bins[16] = {-2.4, -2.1, -1.6, -1.2, -1.04, -0.9, -0.3, -0.2,  0.2, 0.3, 0.9, 1.04, 1.2, 1.6, 2.1, 2.4};
double iso_bins[12] = { 0  , 0.02, 0.04, 0.06, 0.08,  0.1, 0.12, 0.16, 0.2, 0.3, 0.6, 1   };


/// TAG-DEFINITION: 
std::string isofilterTag  = "hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09::HLT";

/// for PROMPT-MUONS   (close-by and far-away) 
std::string L1filter      = "hltL1fL1sMu22or25L1Filtered0::TEST"; 
std::string L2filter      = "hltL2fL1sMu22or25L1f0L2Filtered10Q::TEST";
std::string L3filter      = "hltL3fL1sMu22Or25L1f0L2f10QL3Filtered27Q::TEST"; 

/* 
/// for DIMUONS
std::string L1filter      = "hltL1fL1sDoubleMu155L1Filtered0::TEST"; 
std::string L2filter      = "hltL2pfL1sDoubleMu155L1f0L2PreFiltered0::TEST";
std::string L3filter      = "hltL3fL1DoubleMu155fFiltered17::TEST"; 

/// for DIMUONS (JPSI)
std::string L1filter      = "hltL1fL1sL1sDoubleMu4SQOSdRMax1p2L1Filtered0::TEST"; 
std::string L2filter      = "hltL2fL1sL1sDoubleMu4SQOSdRMax1p2L1f0L2PreFiltered0::TEST";
std::string L3filter      = "hltDimuon25JpsiL3fL3Filtered::TEST"; 

/// for DISPLACED: 
std::string L1filter      = "hltL1fDimuonL1Filtered0::TEST"; 
std::string L2filter      = "hltL2fDimuonL1f0L2NoVtxFiltered16::TEST";
std::string L3filter      = "hltL3fDimuonL1f0L2NVf16L3NoFiltersNoVtxFiltered43::TEST"; 

std::string L1filter      = "hltDimuon3L1Filtered0::TEST"; 
std::string L2filter      = "hltDimuon3L2PreFiltered0::TEST";
std::string L3filter      = "hltDoubleMu3L3FilteredNoVtx::TEST"; 
*/
// ******************************************
//       T&P definitions                    *
//                                          *
std::string thepassfilter  = L3filter;
//std::string theprobefilter = L1filter; 
float offlinePtCut         = 28.;
//                                          *
//                                          *
// ******************************************

void readNtuplesNoFilter_forAN_isMC(TString inputfilename="results.root", int flavor=Sig::Prompt,  std::string effmeasured="IterL3_NOHP_NOL1"){

  bool doingL1 = effmeasured.find("L1OverOffline") != std::string::npos; 
  
  TFile* outfile = TFile::Open(Form("%s_efficiency.root", effmeasured.c_str()),"RECREATE");
  std::cout << "output file: " << outfile -> GetName() << std::endl;

  //Create histograms  
  TH1F* dimuon_mass             = new TH1F("h_dimuon_mass"          ,"dimuon_mass"      , 1500,  0,  150 );
  TH1F* tagiso                  = new TH1F("h_tagiso"               ,"tagiso"           ,  100,  0,  1   );
  TH1F* tagMuonPt               = new TH1F("h_tagMuonPt"            ,"tagMuonPt"        ,  150,  0,  150 );
  TH1F* nvtx_event              = new TH1F("h_nvtx_event"           ,"nvtx_event"       ,   60,  0,   60 );
  
  // L1 over GEN:
  TEfficiency* L1muonPt         = new TEfficiency("L1muonPt"        ,"L1muonPt"         ,   16,  pt_bins ); 
  TEfficiency* L1muonEta        = new TEfficiency("L1muonEta"       ,"L1muonEta"        ,   15, eta_bins );
  TEfficiency* L1muonPhi        = new TEfficiency("L1muonPhi"       ,"L1muonPhi"        ,   20, -3.2, 3.2);
  TEfficiency* L1muonEff        = new TEfficiency("L1muonEff"       ,"L1muonEff"        ,    1,   0., 1.0);
  								    
  // L2 over L1: 						    
  TEfficiency* L2muonPt         = new TEfficiency("L2muonPt"        ,"L2muonPt"         ,   16,  pt_bins ); 
  TEfficiency* L2muonEta        = new TEfficiency("L2muonEta"       ,"L2muonEta"        ,   15, eta_bins );
  TEfficiency* L2muonPhi        = new TEfficiency("L2muonPhi"       ,"L2muonPhi"        ,   20, -3.2, 3.2);
  TEfficiency* L2muonEff        = new TEfficiency("L2muonEff"       ,"L2muonEff"        ,    1,   0., 1.0);

  // L3 over L1:
  TEfficiency* muonPt           = new TEfficiency("muonPt"          ,"muonPt"           ,   16,  pt_bins ); 
  TEfficiency* muonPtTurnOn     = new TEfficiency("muonPtTurnOn"    ,"muonPtTurnOn"     ,   16,  pt_bins ); 
  TEfficiency* muonEta          = new TEfficiency("muonEta"         ,"muonEta"          ,   15, eta_bins );
  TEfficiency* muonPhi          = new TEfficiency("muonPhi"         ,"muonPhi"          ,   20, -3.2, 3.2);
  TEfficiency* muonEff          = new TEfficiency("muonEff"         ,"muonEff"          ,    1,   0., 1.0);
  TEfficiency* muonDeltaR       = new TEfficiency("muonDeltaR"      ,"muonDeltaR"       ,   30,   0., 3.0);
  TEfficiency* muonDeltaPhi     = new TEfficiency("muonDeltaPhi"    ,"muonDeltaPhi"     ,   30,   0., 3.2); 

  TEfficiency* failingMuonPt    = new TEfficiency("failingMuonPt"   ,"failingMuonPt"    ,   16,  pt_bins ); 
  TEfficiency* failingMuonEta   = new TEfficiency("failingMuonEta"  ,"failingMuonEta"   ,   15, eta_bins );
  TEfficiency* failingMuonPhi   = new TEfficiency("failingMuonPhi"  ,"failingMuonPhi"   ,   20, -3.2, 3.2);
  TEfficiency* failingMuonEff   = new TEfficiency("failingMuonEff"  ,"failingMuonEff"   ,   1 ,   0., 1.0);

  TEfficiency* fakeMuonPt    = new TEfficiency("fakeMuonPt"   ,"fakeMuonPt"    ,   16,  pt_bins ); 
  TEfficiency* fakeMuonEta   = new TEfficiency("fakeMuonEta"  ,"fakeMuonEta"   ,   15, eta_bins );
  TEfficiency* fakeMuonPhi   = new TEfficiency("fakeMuonPhi"  ,"fakeMuonPhi"   ,   20, -3.2, 3.2);
  TEfficiency* fakeMuonEff   = new TEfficiency("fakeMuonEff"  ,"fakeMuonEff"   ,   1 ,   0., 1.0);

  TH1F* PassingProbePt          = new TH1F("h_PassingProbePt"       ,"PassingMuonPt"    ,  16,  pt_bins );
  TH1F* PassingProbeEta         = new TH1F("h_PassingProbeEta"      ,"PassingMuonEta"   ,  15, eta_bins );
  TH1F* PassingProbePhi         = new TH1F("h_PassingProbePhi"      ,"PassingMuonPhi"   ,  20, -3.2, 3.2);
  TH1F* PassingProbeMll         = new TH1F("h_PassingProbeMll"      ,"PassingMuonMll"   ,  40, 60., 120.); 
  TH1F* PassingProbeDeltaR      = new TH1F("h_PassingProbeDeltaR"   ,"PassingMuonDeltaR",  30,  0.,  3.0); 

  TH1F* FailingProbePt          = new TH1F("h_FailingProbePt"       ,"FailingMuonPt"    ,  16,  pt_bins );
  TH1F* FailingProbeEta         = new TH1F("h_FailingProbeEta"      ,"FailingMuonEta"   ,  15, eta_bins );
  TH1F* FailingProbePhi         = new TH1F("h_FailingProbePhi"      ,"FailingMuonPhi"   ,  20, -3.2, 3.2);
  TH1F* FailingProbeMll         = new TH1F("h_FailingProbeMll"      ,"FailingMuonMll"   ,  40,  60., 120.);
  TH1F* FailingProbeDeltaR      = new TH1F("h_FailingProbeDeltaR"   ,"FailingMuonDeltaR",  30,  0.,  3.0); 

  // Di-muon efficiencies 
  TEfficiency* diMuonPt         = new TEfficiency("diMuonPt"      ,"diMuonPt"       ,   16,  pt_bins  , 16,  pt_bins ); 
  TEfficiency* diMuonEta        = new TEfficiency("diMuonEta"     ,"diMuonEta"      ,   15, eta_bins  , 15, eta_bins );
  TEfficiency* diMuonPhi        = new TEfficiency("diMuonPhi"     ,"diMuonPhi"      ,   20, -3.2, 3.2 , 20, -3.2, 3.2);
  TEfficiency* diMuonEff        = new TEfficiency("diMuonEff"     ,"diMuonEff"      ,    1,   0., 1.0);
  TEfficiency* diMuonDeltaR     = new TEfficiency("diMuonDeltaR"  ,"diMuonDeltaR"   ,   30,   0., 3.0);
  TEfficiency* diMuonLeadPt     = new TEfficiency("diMuonLeadPt"  ,"diMuonLeadPt"   ,   16,  pt_bins ); 
  TEfficiency* diMuonLeadEta    = new TEfficiency("diMuonLeadEta" ,"diMuonLeadEta"  ,   15,  eta_bins );
  TEfficiency* diMuonLeadPhi    = new TEfficiency("diMuonLeadPhi" ,"diMuonLeadPhi"  ,   20, -3.2, 3.2);
  TEfficiency* diMuonTrailPt    = new TEfficiency("diMuonTrailPt" ,"diMuonTrailPt"  ,   16,  pt_bins ); 
  TEfficiency* diMuonTrailEta   = new TEfficiency("diMuonTrailEta","diMuonTrailEta" ,   15, eta_bins );
  TEfficiency* diMuonTrailPhi   = new TEfficiency("diMuonTrailPhi","diMuonTrailPhi" ,   20, -3.2, 3.2);
  
  TEfficiency* nvtx             = new TEfficiency("nvtx"             ,"nvtx"             ,   60,    0,  60);
  TEfficiency* nvtx_barrel      = new TEfficiency("nvtx_barrel"      ,"nvtx_barrel"      ,   60,    0,  60);
  TEfficiency* nvtx_endcap      = new TEfficiency("nvtx_endcap"      ,"nvtx_endcap"      ,   60,    0,  60);   
 
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

  bool flagfile = false;
  offlinePtCut = getLeadingPtCut(flavor);
  float ptcut1 = getLeadingPtCut(flavor);
  float ptcut2 = getTrailingPtCut(flavor);
      
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)     {
    Int_t IgetEvent   = tree   -> GetEvent(eventNo);
    if (!debug) printProgBar((int)(eventNo*100./nentries));
    if (debug && eventNo==100) break;

    nvtx_event-> Fill( ev -> nVtx   ); 
    unsigned int nmuons = ev->genParticles.size();
    for (int imu = 0; imu < nmuons; imu++){ 
      if (! selectTagMuon(ev->genParticles.at(imu), tagiso)) continue; 
      if (! matchMuon(ev->genParticles.at(imu), ev-> hltTag.objects, isofilterTag)) continue;
      
      tagMuonPt -> Fill ( ev->genParticles.at(imu).pt) ; 
      
      for (int jmu = 0; jmu < nmuons; jmu++){
	bool passL1 = false;
	bool passL2 = false;
	bool passL3 = false;
	
	// Select the probe muon  & match to the probe:  
	if (!selectProbeMuon(ev->genParticles.at(jmu), ev->genParticles.at(imu), dimuon_mass)) continue;
	
	// Check L1 passing: 
	if (matchMuonWithL1(ev->genParticles.at(jmu),ev->L1muons))   passL1=true;
//
//	// Check L2 passing: 
	if (matchMuonWithHLT(ev->genParticles.at(jmu),ev->L2muons))  passL2=true;
//
//	// Check L3 passing: 
	if (matchMuonWithHLT(ev->genParticles.at(jmu),ev->hltmuons)) passL3=true;
	
	L1muonPt     -> Fill( passL1          , ev->genParticles.at(jmu).pt ); 
	L2muonPt     -> Fill( passL2 && passL1, ev->genParticles.at(jmu).pt ); 
	muonPtTurnOn -> Fill( passL3 && passL1, ev->genParticles.at(jmu).pt ); 
	
	if (ev->genParticles.at(jmu).pt < ptcut1) continue;
	L1muonEta    -> Fill( passL1, ev->genParticles.at(jmu).eta);
	L1muonPhi    -> Fill( passL1, ev->genParticles.at(jmu).phi);
	L1muonEff    -> Fill( passL1, 0.5                    );
	
	if (!passL1) continue;  // exit the loop if there is no L1
	
	// Fill L2:
	L2muonEta    -> Fill( passL2, ev->genParticles.at(jmu).eta);
	L2muonPhi    -> Fill( passL2, ev->genParticles.at(jmu).phi);
	L2muonEff    -> Fill( passL2, 0.5                    );

	// Fill L3:
	TLorentzVector mu1, mu2;
	mu1.SetPtEtaPhiM (ev->genParticles.at(imu).pt,ev->genParticles.at(imu).eta,ev->genParticles.at(imu).phi, muonmass); 
	mu2.SetPtEtaPhiM (ev->genParticles.at(jmu).pt,ev->genParticles.at(jmu).eta,ev->genParticles.at(jmu).phi, muonmass);
	double mumumass = (mu1 + mu2).M();
	double DeltaR   = mu1.DeltaR(mu2);
	double DeltaPhi = mu1.DeltaPhi(mu2);

	if (passL3) { 
	  PassingProbePt    -> Fill( ev->genParticles.at(jmu).pt  );
	  PassingProbeEta   -> Fill( ev->genParticles.at(jmu).eta );
	  PassingProbePhi   -> Fill( ev->genParticles.at(jmu).phi );
	  PassingProbeMll   -> Fill( mumumass                );
	  PassingProbeDeltaR-> Fill( DeltaR );
	}	      
	else {       
	  FailingProbePt    -> Fill( ev->genParticles.at(jmu).pt  );
	  FailingProbeEta   -> Fill( ev->genParticles.at(jmu).eta );
	  FailingProbePhi   -> Fill( ev->genParticles.at(jmu).phi );
	  FailingProbeMll   -> Fill( mumumass                );
	  FailingProbeDeltaR-> Fill( DeltaR );
	}
	
	muonPt       -> Fill( passL3, ev->genParticles.at(jmu).pt );
	muonEta      -> Fill( passL3, ev->genParticles.at(jmu).eta);
	muonPhi      -> Fill( passL3, ev->genParticles.at(jmu).phi);
	muonEff      -> Fill( passL3, 0.5                    );
	muonDeltaR   -> Fill( passL3, DeltaR);
	muonDeltaPhi -> Fill( passL3, DeltaPhi);
	
	failingMuonPt  -> Fill( !passL3, ev->genParticles.at(jmu).pt );
	failingMuonEta -> Fill( !passL3, ev->genParticles.at(jmu).eta);
	failingMuonPhi -> Fill( !passL3, ev->genParticles.at(jmu).phi);
	failingMuonEff -> Fill( !passL3, 0.5                    );
      } // nmuons
    }
    
    /// Now FOR DIMUONS (Using HLT info):
    // both should pass the L1 filter! 
    continue;
    if (debug) cout << "fired L1?" << endl;
    if (ev->L1muons.size()>1) { 
      if (debug) cout << "YES!" << endl;      
      for (int imu=0; imu<ev->hltmuons.size(); imu++){ 
	bool fake = false;
	if (!matchHLTMuonWithGen(ev->hltmuons.at(imu), ev->genParticles)) fake = true;
	fakeMuonPt  -> Fill( fake, ev->hltmuons.at(imu).pt );
	if (ev->hltmuons.at(imu).pt > ptcut1) { 
	  fakeMuonEta -> Fill( fake, ev->hltmuons.at(imu).eta);
	  fakeMuonPhi -> Fill( fake, ev->hltmuons.at(imu).phi);
	  fakeMuonEff -> Fill( fake, 0.5                     );
	}
	for (int jmu=imu+1; jmu<ev->hltmuons.size(); jmu++){ 
	  bool pass = false;
	  if (matchHLTMuonWithGen(ev->hltmuons.at(imu),ev->genParticles) && 
	      matchHLTMuonWithGen(ev->hltmuons.at(jmu),ev->genParticles)) pass = true;
	  
	  double DeltaR = deltaR(ev->hltmuons.at(imu).eta,ev->hltmuons.at(imu).phi,ev->hltmuons.at(jmu).eta,ev->hltmuons.at(jmu).phi);
	  diMuonPt       -> Fill( pass, ev->hltmuons.at(imu).pt , ev->hltmuons.at(jmu).pt  ); 
	  diMuonLeadPt   -> Fill( pass, ev->hltmuons.at(imu).pt ); 
	  diMuonTrailPt  -> Fill( pass, ev->hltmuons.at(jmu).pt ); 
	  
	  if (ev->hltmuons.at(imu).pt < ptcut1) continue;
	  if (ev->hltmuons.at(jmu).pt < ptcut2) continue;
	  diMuonEta      -> Fill( pass, ev->hltmuons.at(imu).eta, ev->hltmuons.at(jmu).eta );
	  diMuonPhi      -> Fill( pass, ev->hltmuons.at(imu).phi, ev->hltmuons.at(jmu).phi );
	  diMuonDeltaR   -> Fill( pass, DeltaR);
	  diMuonEff      -> Fill( pass, 0.5);
	  
	  diMuonLeadEta  -> Fill( pass, ev->hltmuons.at(imu).eta ); 
	  diMuonTrailEta -> Fill( pass, ev->hltmuons.at(jmu).eta ); 
	  diMuonLeadPhi  -> Fill( pass, ev->hltmuons.at(imu).phi ); 
	  diMuonTrailPhi -> Fill( pass, ev->hltmuons.at(jmu).phi ); 
	}
      }
      /*
	for (int imu = 0; imu < nmuons; imu++){
	if (!selectMuon(ev->genParticles.at(imu))) continue; 
	for (int jmu = imu+1; jmu < nmuons; jmu++){
	if (!selectMuon(ev->genParticles.at(jmu))) continue; 
	bool pass = false;
	if (matchMuonWithL3(ev->genParticles.at(jmu),ev->hltmuons) && 
	matchMuonWithL3(ev->genParticles.at(imu),ev->hltmuons)) pass = true;
	
	if (ev->genParticles.at(imu).pt  == ev->genParticles.at(jmu).pt)  continue;
	if (ev->genParticles.at(imu).eta == ev->genParticles.at(jmu).eta) continue;
	if (ev->genParticles.at(imu).phi == ev->genParticles.at(jmu).phi) continue;
	
	double DeltaR   = deltaR(ev->genParticles.at(imu).eta,ev->genParticles.at(imu).phi,ev->genParticles.at(jmu).eta,ev->genParticles.at(jmu).phi);
	h_diMuonDeltaR -> Fill( DeltaR );
	
	diMuonPt       -> Fill( pass, ev->genParticles.at(imu).pt ); 
	diMuonEta      -> Fill( pass, ev->genParticles.at(imu).eta);
	diMuonPhi      -> Fill( pass, ev->genParticles.at(imu).phi);
	diMuonPt       -> Fill( pass, ev->genParticles.at(jmu).pt );
	diMuonEta      -> Fill( pass, ev->genParticles.at(jmu).eta);
	diMuonPhi      -> Fill( pass, ev->genParticles.at(jmu).phi);
	diMuonDeltaR   -> Fill( pass, DeltaR);
	diMuonDeltaR   -> Fill( pass, DeltaR);
	diMuonEff      -> Fill( pass, 0.5);
	diMuonEff      -> Fill( pass, 0.5);
	}
	}
      */
    }
  }  
  
  //Writing the histograms in a file.
  outfile           -> cd();
  tagMuonPt         -> Write();

  L1muonPt            -> Write();
  L1muonEta           -> Write();
  L1muonPhi           -> Write();
  L1muonEff           -> Write();
  
  L2muonPt            -> Write();
  L2muonEta           -> Write();
  L2muonPhi           -> Write();
  L2muonEff           -> Write();
  
  muonPtTurnOn      -> Write();
  muonPt            -> Write();
  muonEta           -> Write();
  muonPhi           -> Write();
  muonEff           -> Write();
  muonDeltaR        -> Write();
  muonDeltaPhi      -> Write();

  failingMuonPt   -> Write();
  failingMuonEta  -> Write();
  failingMuonPhi  -> Write();
  failingMuonEff  -> Write();

  fakeMuonPt   -> Write();
  fakeMuonEta  -> Write();
  fakeMuonPhi  -> Write();
  fakeMuonEff  -> Write();

  PassingProbePt  -> Write();
  PassingProbeEta -> Write();
  PassingProbePhi -> Write();
  PassingProbeMll -> Write();
  PassingProbeDeltaR -> Write();
  
  FailingProbePt  -> Write();
  FailingProbeEta -> Write();
  FailingProbePhi -> Write();
  FailingProbeMll -> Write();
  FailingProbeDeltaR -> Write();

  /// Per-event (di-muon) efficiency.
  diMuonPt         -> Write();
  diMuonEta        -> Write();
  diMuonPhi        -> Write();
  diMuonEff        -> Write();
  diMuonDeltaR     -> Write();
  diMuonLeadPt     -> Write(); 
  diMuonLeadEta    -> Write();
  diMuonLeadPhi    -> Write();
  diMuonTrailPt    -> Write(); 
  diMuonTrailEta   -> Write();
  diMuonTrailPhi   -> Write();

  nvtx_event       -> Write();
  nvtx             -> Write();
  
  dimuon_mass      -> Write();
  tagiso           -> Write();
  
  outfile          -> Close();  
  
  return;
}
bool firedL1( std::vector<HLTObjCand> toc, std::string L1FilterName){ 
  int ntoc = toc.size();
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) { 
    if (debug) cout << it->filterTag << " = " << L1FilterName << endl;
    if ( it->filterTag.compare(L1FilterName) == 0) return true;
  }
  return false;
}
bool matchMuon(GenParticleCand mu, std::vector<HLTObjCand> toc, std::string tagFilterName){

  bool match = false;
  int ntoc = toc.size();

  float minDR = 0.1; 
  if (tagFilterName.find("L1fL1") != std::string::npos) minDR = 1.0;
  float theDR = 100;
  for ( std::vector<HLTObjCand>::const_iterator it = toc.begin(); it != toc.end(); ++it ) { 
    //if (debug) cout << tagFilterName << " = " << it->filterTag << " ? " << endl;
    if ( it->filterTag.compare(tagFilterName) == 0) { 
      //      if (debug) cout << " YES" << endl;
      theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi);
      if (theDR < minDR){
        minDR = theDR;
        match = true;
      }
    }
  }
  
  return match;
}

bool selectTagMuon(GenParticleCand mu, TH1F* tagh){
  if (!( fabs(mu.pdgId) == 13)) return false;
  if (!( mu.pt         > offlinePtCut)) return false; 
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  
  //add isolation cut
  tagh -> Fill(0.);
  return true;
}
float getLeadingPtCut(int signature){ 
  float ptcut = 0.;
  if (signature == Sig::Prompt) ptcut = 29.;
  if (signature == Sig::DiMuon) ptcut = 18.;
  if (signature == Sig::LowPt ) ptcut = 0.;
  return ptcut;
}

float getTrailingPtCut(int signature){ 
  float ptcut = 0.;
  if (signature == Sig::Prompt) ptcut = 27.;
  if (signature == Sig::DiMuon) ptcut = 8. ;
  if (signature == Sig::LowPt ) ptcut = 0. ;
  return ptcut;
}

bool selectMuon(GenParticleCand mu){  
  if (!( fabs(mu.pdgId) == 13))           return false;
  if (!( mu.pt         > offlinePtCut  )) return false; 
  if (!( fabs(mu.eta)  < 2.4 ))           return false;
  return true;
}


//select the probe muon
bool selectProbeMuon(GenParticleCand mu, GenParticleCand tagMu, TH1F* dimuon_mass){
  
  if (mu.pt  == tagMu.pt  && 
      mu.eta == tagMu.eta &&
      mu.phi == tagMu.phi ) 
    return false;
  
  if (!( fabs(mu.eta)  < 2.4 )) return false; 
  if (mu.pdgId * tagMu.pdgId > 0) return false;
  
  TLorentzVector mu1, mu2;
  mu1.SetPtEtaPhiM (mu.pt   , mu.eta   , mu.phi   , muonmass);
  mu2.SetPtEtaPhiM (tagMu.pt, tagMu.eta, tagMu.phi, muonmass);
  double mumumass = (mu1 + mu2).M();
  dimuon_mass -> Fill(mumumass); 
  if (! (mumumass > 80. && mumumass < 110. )) return false;
  
  return true;
}

bool matchMuonWithHLT(GenParticleCand mu, std::vector<HLTMuonCand> L3cands){

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
bool matchMuonWithL1 (GenParticleCand mu, std::vector<L1MuonCand> L1cands){

  bool match = false;
  float minDR = 1.0;
  float theDR = 100;
  for ( std::vector<L1MuonCand>::const_iterator it = L1cands.begin(); it != L1cands.end(); ++it ) { 
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi); 
    if (theDR < minDR){ 
      minDR = theDR;
      match = true;
    }
  }
  return match;
}

bool matchHLTMuonWithGen(HLTMuonCand mu, std::vector<GenParticleCand> Gencands){

  bool match = false;
  float minDR = 0.1;
  float theDR = 100;
  for ( std::vector<GenParticleCand>::const_iterator it = Gencands.begin(); it != Gencands.end(); ++it ) { 
    theDR = deltaR(it -> eta, it -> phi, mu.eta, mu.phi); 
    if (theDR < minDR){ 
      minDR = theDR;
      match = true;
    }
  }
  return match;
}

bool matchL1MuonWithGen(L1MuonCand mu, std::vector<GenParticleCand> Gencands){

  bool match = false;
  float minDR = 1.0;
  float theDR = 100;
  for ( std::vector<GenParticleCand>::const_iterator it = Gencands.begin(); it != Gencands.end(); ++it ) { 
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
