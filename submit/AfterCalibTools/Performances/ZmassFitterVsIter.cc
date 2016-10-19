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
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

#include <fstream>
#include <iostream>        

using namespace RooFit;
using namespace std;  

// ============================================
// to be modified:
static const Int_t NLOOP = 10;  
Int_t MINmass= 70;
Int_t MAXmass= 110;
// ============================================


void doAll() {
  
  Int_t nloop = NLOOP;

  // TFile
  cout << endl; cout << "File and Tree" << endl;
  TFile *tfile = TFile::Open("allEpsilonsRereco.root");

  // Tree
  TTree *ttree = (TTree*)tfile->Get("myTree"); 
  float zMass;
  int iteration;
  int isEBEB;
  int isEEEE;
  ttree->SetBranchAddress("zMass",&zMass);
  ttree->SetBranchAddress("iteration",&iteration);
  ttree->SetBranchAddress("isEBEB",&isEBEB);
  ttree->SetBranchAddress("isEEEE",&isEEEE);

  // RooRealVar
  cout <<  endl; cout << "RooRealVars" << endl;
  RooRealVar massRRV("zmassRRV", "m(e^{+}e^{-})", MINmass, MAXmass, "GeV");
  RooArgSet zMassArgSet(massRRV);

  // to save fit results  
  vector<float> v_iter, v_mean, v_sigma;
  vector<float> v_iterE, v_meanE, v_sigmaE;

  // Loop over iteration and fit
  cout << endl; cout << "preparing dataset" << endl;

  for (int iLoop=0; iLoop<nloop; ++iLoop) {

    RooDataSet* dataR = new RooDataSet("dataR", "ntuple parameters", zMassArgSet);

    for (int i = 0; i < ttree->GetEntries(); i++) {
      if(i%5000000==0) cout << "Loop " << iLoop << ", Processing Event " << i << endl;
      ttree->GetEntry(i);

      // reduction - chiara
      if (isEEEE!=1) continue;
      if (iteration!=iLoop) continue;

      // mass selection
      float zMassRRV = zMass;
      if (zMassRRV<MINmass || zMassRRV>MAXmass) continue;

      zMassArgSet.setRealValue("zmassRRV", zMassRRV);

      dataR->add(zMassArgSet);
    }
    
    cout << endl; cout << "reduced dataset size for this loop: " << dataR->numEntries() << endl;
    
    // Binned dataset
    massRRV.setBins(40);   
    RooDataHist *bdataR = new RooDataHist("data_binned","data_binned", zMassArgSet, *dataR);    
    
    // Fit: dCB
    RooRealVar CBmean ("CBmean", "CBmean", 0.,-3.,1.);             // scala ok
    // RooRealVar CBmean ("CBmean", "CBmean", 0.,-3.,15.);            // scala +10%
    // RooRealVar CBmean ("CBmean", "CBmean", 0.,-20.,1.);            // scala -10%
    //RooRealVar CBsigma("CBsigma","CBsigma", 0.8, 0.01, 3.);        // no miscalib EBEB
    RooRealVar CBsigma("CBsigma","CBsigma", 0.8, 0.01, 10.);        // no miscalib EEEE
    // RooRealVar CBsigma("CBsigma","CBsigma", 0.8, 0.01, 7.);
    RooRealVar CBalpha1("CBalpha1","CBalpha1",1.,0.,2.);
    RooRealVar CBn1("CBn1","CBn1",1.5,0.,10.);
    RooRealVar CBalpha2("CBalpha2","CBalpha2",1.,0.,3.);
    RooRealVar CBn2("CBn2","CBn2",1.5,0.,10.);
    RooDoubleCB ResCB("ResCB","ResCB", massRRV, CBmean, CBsigma, CBalpha1, CBn1, CBalpha2, CBn2);
    
    // Fit: BW 
    RooRealVar meanBW("massBW","massBW",91.1876);
    RooRealVar sigmaBW("widthBW","widthBW",2.4952);
    RooBreitWigner zBW("BW","BW", massRRV, meanBW, sigmaBW);
    
    // Fit: Convolution
    RooFFTConvPdf* ConvolutedRes = new RooFFTConvPdf("conv","conv", massRRV, zBW, ResCB);
      
    // Now fitting
    RooFitResult* fitres = ConvolutedRes->fitTo(*bdataR, SumW2Error(kFALSE), Range(MINmass,MAXmass), RooFit::Save(kTRUE));
    fitres->SetName("fitres");  
    fitres->Print("V");
    
    // Saving the fit outcome
    v_iter.push_back(iLoop);
    v_iterE.push_back(0);
    v_mean.push_back(CBmean.getVal()+91.1876);
    v_meanE.push_back(CBmean.getError());
    v_sigma.push_back(CBsigma.getVal());
    v_sigmaE.push_back(CBsigma.getError());
    
    // Now plotting
    cout << endl; cout << "Now plotting" << endl;
    TCanvas* c = new TCanvas("c","Invariant Mass Fit", 0,0,800,600);    
    RooPlot* plot = massRRV.frame();     
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

    c->SaveAs(TString::Format("zmass_iter%d.png",iLoop));
    c->SaveAs(TString::Format("zmass_iter%d.pdf",iLoop));

    delete tex;
    delete c;
    delete ConvolutedRes;
    delete bdataR;

  }  // loop over iterations

  // graphs with trend
  TGraphErrors *gMean  = new TGraphErrors(NLOOP, &v_iter[0], &v_mean[0],  &v_iterE[0], &v_meanE[0]);
  TGraphErrors *gSigma = new TGraphErrors(NLOOP, &v_iter[0], &v_sigma[0], &v_iterE[0], &v_sigmaE[0]);
  TFile fileOut("outTrends.root","RECREATE");
  gMean->Write("gMean");
  gSigma->Write("gSigma");
  fileOut.Close();
  
  delete gMean;
  delete gSigma;
}

void runfits() {

  doAll();   

  return;
}

