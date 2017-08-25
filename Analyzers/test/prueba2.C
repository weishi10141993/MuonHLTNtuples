#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
#include "TGaxis.h"

void prueba2(){
   
  TCanvas *c1 = new TCanvas("myCanvas","myCanvas",0,600);
  TCanvas *c2 = new TCanvas("myCanvas1","myCanvas1",0,600);
  TCanvas *c3 = new TCanvas("myCanvas2","myCanvas2",0,600);


  TFile *f1 = TFile::Open("Cascade_Mu_HltTrack_Match.root", "READ");
  TFile *f2 = TFile::Open("SeedFix_Mu_HltTrack_Match.root", "READ");
  TFile *f3 = TFile::Open("Navigation_Mu_HltTrack_Match.root", "READ");


  TH1F *h1  = (TH1F*) f1 -> Get("h_PtHLT_barrel");
  TH1F *h2  = (TH1F*) f2 -> Get("h_PtHLT_barrel");
  TH1F *h3  = (TH1F*) f3 -> Get("h_PtHLT_barrel");

  TH1F *h4  = (TH1F*) f1 -> Get("h_PtHLT_int");
  TH1F *h5  = (TH1F*) f2 -> Get("h_PtHLT_int");
  TH1F *h6  = (TH1F*) f3 -> Get("h_PtHLT_int");

  TH1F *h7  = (TH1F*) f1 -> Get("h_PtHLT_endcap");
  TH1F *h8  = (TH1F*) f2 -> Get("h_PtHLT_endcap");
  TH1F *h9  = (TH1F*) f3 -> Get("h_PtHLT_endcap");

  auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
  legend1->AddEntry(h1,"Cascade+TkMu","l");
  legend1->AddEntry(h2,"HitBased+Grouped","l");
  legend1->AddEntry(h3,"Hitless+HitBased","l");

  int nl1 = h1 -> GetEntries();
  int nl2 = h2 -> GetEntries();
  int nl3 = h3 -> GetEntries();
  int nl4 = h4 -> GetEntries();
  int nl5 = h5 -> GetEntries();
  int nl6 = h6 -> GetEntries();
  int nl7 = h7 -> GetEntries();
  int nl8 = h8-> GetEntries();
  int nl9 = h9-> GetEntries();

  h1 -> Scale(1.0/nl1);
  h2 -> Scale(1.0/nl2);
  h3 -> Scale(1.0/nl3);
  h4 -> Scale(1.0/nl4);
  h5 -> Scale(1.0/nl5);
  h6 -> Scale(1.0/nl6);
  h7 -> Scale(1.0/nl7);
  h8 -> Scale(1.0/nl8);
  h9 -> Scale(1.0/nl9);

  //c1 -> Divide(2,1);
  //h2 -> SetLineWidth (2);
  //h2 -> SetAxisRange( -0.5, 6.5, "X");

  c1 -> cd();
  h3 -> Draw("HIST");
  h3 -> SetLineColor(kGreen);
  h3 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h3 -> SetAxisRange( -0.15, 0.15, "X");

  h1 -> Draw("HIST same");
  h1 -> SetLineColor(kBlack);
  h1 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h1 -> SetAxisRange( -0.15, 0.15, "X");

  h2 -> Draw("HIST same");
  h2 -> SetLineColor(kRed);
  h2 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h2 -> SetAxisRange( -0.15, 0.15, "X");
  legend1 -> Draw();
  c1 -> Update();
  c1 -> Draw("");




  c2 -> cd();
  h6 -> Draw("HIST");
  h6 -> SetLineColor(kGreen);
  h6 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h6 -> SetAxisRange( -0.15, 0.15, "X");

  h5 -> Draw("HIST same");
  h5 -> SetLineColor(kRed);
  h5 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h5 -> SetAxisRange( -0.15, 0.15, "X");

  h4 -> Draw("HIST same");
  h4 -> SetLineColor(kBlack);
  h4 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h4 -> SetAxisRange( -0.15, 0.15, "X");
  legend1 -> Draw();
  c2 -> Update();
  c2 -> Draw("");


  c3 -> cd();
  h8 -> Draw("HIST");
  h8 -> SetLineColor(kRed);
  h8 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h8 -> SetAxisRange( -0.15, 0.15, "X");

  h9 -> Draw("HIST same");
  h9 -> SetLineColor(kGreen);
  h9 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h9 -> SetAxisRange( -0.15, 0.15, "X");;

  h7 -> Draw("HIST Same");
  h7 -> SetLineColor(kBlack);
  h7 -> GetXaxis()->SetTitle("(pT muon - pT HLT muon)/(pT muon)");
  h7 -> SetAxisRange( -0.15, 0.15, "X");
  legend1 -> Draw();
  c3 -> Update();
  c3 -> Draw("");

}
