#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
#include "TGaxis.h"

int prueba(){
   
  TCanvas *c1 = new TCanvas("myCanvas","myCanvas",600,600);
  TFile* inputfile = TFile::Open("L3_efficiency_newIter_old__pre.root", "READ");
  TFile* outfile   = TFile::Open("prueba.root", "RECREATE");

  TEfficiency* global  = (TEfficiency*) inputfile -> Get("muonPixHit");
  TEfficiency* barrel  = (TEfficiency*) inputfile -> Get("muonPixHit_barrel");
  TEfficiency* overlap = (TEfficiency*) inputfile -> Get("muonPixHit_int");
  TEfficiency* endcap  = (TEfficiency*) inputfile -> Get("muonPixHit_endcap");


  c1 -> cd(1);

  global ->Draw();
  global ->SetLineColor(kBlack);
  barrel ->Draw("same");
  barrel ->SetLineColor(kBlue);
  overlap->Draw("same");
  overlap->SetLineColor(kRed);
  endcap ->Draw("same");
  endcap ->SetLineColor(kGreen);

  c1 -> Update();


  return;
}
