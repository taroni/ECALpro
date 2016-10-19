#include "TStyle.h"
#include "TFile.h"
#include <TGraphErrors.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2F.h"
#include <iostream>

using namespace std;

void compareTrends() {

  // taking inputs                                                                                                                           
  TFile *infile[2];
  infile[0] = new TFile("outTrendsPromptEBEB.root");
  infile[1] = new TFile("outTrendsRerecoEBEB.root");
  for (int ii=0; ii<2; ii++) {
    if (!infile[ii]) cout << "File " << infile[ii] << " not existing" << endl;
  }

  TGraphErrors *Gmean[2], *Gsigma[2];   
  for (int ii=0; ii<2; ii++) {     
    Gmean[ii]  = (TGraphErrors*)infile[ii]->Get(TString::Format("gMean"));
    Gsigma[ii] = (TGraphErrors*)infile[ii]->Get(TString::Format("gSigma"));
    Gmean[ii]  -> SetMarkerStyle(20);    
    Gsigma[ii] -> SetMarkerStyle(20);
    Gmean[ii]  -> SetMarkerColor(ii+1);    
    Gsigma[ii] -> SetMarkerColor(ii+1);
  }
  

  // ---------------------------------------------------------------------
  TLegend* leg = new TLegend(0.15, 0.10, 0.75, 0.25); 
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(Gmean[0], "prompt reco", "p");  
  leg->AddEntry(Gmean[1], "re-reco", "p");  

  gStyle->SetOptStat(0); 

  TH2F *myH_mean = new TH2F("myH_mean","",100,-0.5,13.5,100,90.,93.);
  myH_mean->SetTitle("");
  myH_mean->GetXaxis()->SetTitle("iteration");
  myH_mean->GetYaxis()->SetTitle("CB mean [GeV]");

  TCanvas c1m("c1m","cat0",1);
  myH_mean -> Draw();
  for (int ii=0; ii<2; ii++) Gmean[ii] ->Draw("Psame");
  leg->Draw();  
  c1m.SaveAs("mass.png");

  // ---------------------------------------------------------------------
  TH2F *myH_sigma = new TH2F("myH_sigma","myH_sigma",100,-0.5,13.5,100,1.3,2.3);
  //TH2F *myH_sigma = new TH2F("myH_sigma","myH_sigma",100,-0.5,13.5,100,2.5,3.5);
  myH_sigma->SetTitle("");
  myH_sigma->GetXaxis()->SetTitle("iteration");
  myH_sigma->GetYaxis()->SetTitle("CB sigma [GeV]");

  TCanvas c1s("c1s","cat0",1);
  myH_sigma -> Draw();
  for (int ii=0; ii<2; ii++) Gsigma[ii] ->Draw("Psame");    
  leg->Draw();  
  c1s.SaveAs("sigma.png");
}
