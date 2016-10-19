#include "RooRealVar.h"
#include "RooBinning.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

#include <fstream>
#include <iostream>

using namespace RooFit;
using namespace std;

// ============================================
// to be modified:
Int_t MINmass= 70;
Int_t MAXmass= 110;
// ============================================

void AddSigData() {

  // -----------------------------------------------------------------
  
  // Variables
  cout <<  endl; cout << "variables" << endl;
  RooRealVar zMass("zMass", "m(e^{+}e^{-})", MINmass, MAXmass, "GeV");
  RooRealVar isEBEB("isEBEB",    "isEBEB",   -0.1, 1.1,   "");
  RooRealVar isEEEE("isEEEE",    "isEEEE",   -0.1, 1.1,   "");
  RooRealVar iteration("iteration", "iteration", 0, 30, "");
  
  // Dataset from tree
  cout << endl; cout << "dataset" << endl;
  TFile *tfile = TFile::Open("TestZ_epsilonPlots7.root");
  TTree *ttree = (TTree*)tfile->Get("myTree"); 
  RooArgSet zMassArgSet(zMass,isEBEB,isEEEE,iteration);
  RooDataSet* data= new RooDataSet("data","dataset",ttree,zMassArgSet,"");   
  cout << endl; cout << "dataset size: " << data->numEntries() << endl;  

  // Reduced dataset
  cout << endl; cout << "reduction" << endl;
  RooDataSet* dataR = (RooDataSet*) data->reduce(zMass,TString::Format("isEBEB"));
  cout << endl; cout << "reduced dataset size: " << dataR->numEntries() << endl;

  // Binned dataset
  zMass.setBins(40);   
  RooDataHist *bdataR = new RooDataHist("data_binned","data_binned", zMassArgSet, *dataR);    
  


  // -----------------------------------------------------------------
      
  // Fit
  // dCB
  RooRealVar CBmean ("CBmean", "CBmean", 0.,-3.,1.);             // scala ok
  // RooRealVar CBmean ("CBmean", "CBmean", 0.,-3.,15.);            // scala +10%
  // RooRealVar CBmean ("CBmean", "CBmean", 0.,-20.,1.);            // scala -10%
  RooRealVar CBsigma("CBsigma","CBsigma", 0.8, 0.01, 3.);        // no miscalib
  // RooRealVar CBsigma("CBsigma","CBsigma", 0.8, 0.01, 7.);
  RooRealVar CBalpha1("CBalpha1","CBalpha1",1.,0.,2.);
  RooRealVar CBn1("CBn1","CBn1",1.5,0.,10.);
  RooRealVar CBalpha2("CBalpha2","CBalpha2",1.,0.,3.);
  RooRealVar CBn2("CBn2","CBn2",1.5,0.,10.);
  RooDoubleCB ResCB("ResCB","ResCB", zMass, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2);
  
  // BW 
  RooRealVar meanBW("massBW","massBW",91.1876);
  RooRealVar sigmaBW("widthBW","widthBW",2.4952);
  RooBreitWigner zBW("BW","BW", zMass, meanBW, sigmaBW);
  
  // Convolution
  RooFFTConvPdf* ConvolutedRes;
  ConvolutedRes = new RooFFTConvPdf("conv","conv", zMass, zBW, ResCB);
      
  //////////////////////////////////////////////////////////////

  RooFitResult* fitres = ConvolutedRes->fitTo(*bdataR, SumW2Error(kFALSE), Range(MINmass,MAXmass), RooFit::Save(kTRUE));
  fitres->SetName("fitres");  

  cout << endl; cout << "Now plotting" << endl;

  TCanvas* c = new TCanvas("c","Invariant Mass Fit", 0,0,800,600);    
  RooPlot* plot = zMass.frame();     
  dataR->plotOn(plot);  
  ConvolutedRes->plotOn(plot,LineColor(kOrange));  
  plot->Draw();   
  
  // Print Fit Values     
  TLatex *tex = new TLatex();   
  tex->SetNDC();      
  tex->SetTextSize(.1); 
  tex->SetTextFont(132);    
  tex->SetTextSize(0.057);
  tex->DrawLatex(0.65, 0.75, "Z #rightarrow e^{+}e^{-}");    
  tex->SetTextSize(0.030);      
  tex->DrawLatex(0.645, 0.65, Form("CB #DeltaMean = %.3f GeV", CBmean.getVal()));    
  tex->DrawLatex(0.645, 0.60, Form("CB #sigma = %.3f GeV", CBsigma.getVal()));   
  c->Update();                     
  c->SaveAs("zmass.pdf");
  c->SaveAs("zmass.png");
}

void runfits() {

  cout << endl; 
  cout << "Now add data" << endl;
  AddSigData();   

  return;
}

