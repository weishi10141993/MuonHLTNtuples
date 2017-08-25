#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
#include "TGaxis.h"

void plots_same_canvas(){
   
  TCanvas *c1 = new TCanvas("myCanvas","myCanvas",0,600);
  TFile *file = TFile::Open("Cascade_IterL3_Match.root", "READ");

  // TH1F *h1  = (TH1F*) file -> Get("h_hltOI_diff_endcap");
  TH1F *h2  = (TH1F*) file -> Get("h_Pt_delta");
  TH1F *h3  = (TH1F*) file -> Get("h_eta_delta");  
  TH1F *h4  = (TH1F*) file -> Get("h_phi_delta"); 
  TH1F *h5  = (TH1F*) file -> Get("h_chi2_delta");
  TH1F *h6  = (TH1F*) file -> Get("h_dxy_delta");
   TH1F *h7  = (TH1F*) file -> Get("h_dz_delta");

  // TH1F *h8  = (TH1F*) file -> Get("h_FracValTrackHit_delta");
  // TH1F *h9  = (TH1F*) file -> Get("h_PixHit_delta");
  // TH1F *h10  = (TH1F*) file -> Get("h_LayHit_delta");
  // TH1F *h11  = (TH1F*) file -> Get("h_PixLay_delta");

  c1 -> Divide(3,2);
  // c1 -> cd(1);
  // h1 -> Draw();
  c1 -> cd(1);
  h2 -> Draw();
  c1 -> cd(2);
  h3 -> Draw();
  c1 -> cd(3);
  h4 -> Draw();
  c1 -> cd(4);
  h5 -> Draw();
  c1 -> cd(5);
  h6 -> Draw();
   c1 -> cd(6);
  h7 -> Draw();

  // c1 -> cd(1);
  // h8 -> Draw();
  // c1 -> cd(2);
  // h9 -> Draw();
  // c1 -> cd(3);
  // h10 -> Draw();
  // c1 -> cd(4);
  // h11 -> Draw();


  c1 -> Update();

  c1 -> Draw("");

}
