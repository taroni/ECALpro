// system include files
#include <memory>
#include <iostream>
#include <string>

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
#include "TMath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooDataHist.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"

#include "CalibCode/FitEpsilonPlot/interface/FitEpsilonPlotForZ.h"

using std::cout;
using std::endl;
using std::string;

#define MZ 91.1876

using namespace RooFit;

FitEpsilonPlotForZ::FitEpsilonPlotForZ(const edm::ParameterSet& iConfig) {
  
  // now do what ever initialization is needed
  currentIteration_ =  iConfig.getUntrackedParameter<int>("CurrentIteration");
  epsilonPlotFileName_ = iConfig.getUntrackedParameter<std::string>("EpsilonPlotFileName");
  outputDir_ = iConfig.getUntrackedParameter<std::string>("OutputDir");
  outfilename_          = iConfig.getUntrackedParameter<std::string>("OutputFile");
  calibMapPath_ = iConfig.getUntrackedParameter<std::string>("calibMapPath");
  inRangeFit_ = iConfig.getUntrackedParameter<int>("NInFit");
  finRangeFit_ = iConfig.getUntrackedParameter<int>("NFinFit");    
  EEoEB_ = iConfig.getUntrackedParameter<std::string>("EEorEB");
  StoreForTest_ = iConfig.getUntrackedParameter<bool>("StoreForTest","false");
  Barrel_orEndcap_ = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");
  useMassInsteadOfEpsilon_ = iConfig.getUntrackedParameter<bool>("useMassInsteadOfEpsilon",false);

  // setting calibration type
  regionalCalibration_ = &xtalCalib; 
  cout << "FIT_EPSILON_FOR_Z: crosscheck: selected type: " << regionalCalibration_->printType() << endl;

  // retrieving calibration coefficients of the previous iteration
  char fileName[200];
  if(currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
  else if(currentIteration_ > 0) {
    sprintf(fileName,"%s", calibMapPath_.c_str());
    regionalCalibration_->getCalibMap()->loadCalibMapFromFile(fileName);
  }

  // # of regions
  nRegionsEB_ = regionalCalibration_->getCalibMap()->getNRegionsEB() ;
  nRegionsEE_ = regionalCalibration_->getCalibMap()->getNRegionsEE() ;

  // load epsilon from current iter 
  weightedRescaleFactorEB = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEB()];
  weightedRescaleFactorEE = new TH1F*[regionalCalibration_->getCalibMap()->getNRegionsEE()];
  sprintf(fileName,"%s", epsilonPlotFileName_.c_str());
  cout << "FIT_EPSILON: FitEpsilonPlotForZ:: loading epsilon plots from file: " << epsilonPlotFileName_ << endl;
  loadEpsilonPlot(fileName);
}

FitEpsilonPlotForZ::~FitEpsilonPlotForZ() {
  
  delete weightedRescaleFactorEB;
  delete weightedRescaleFactorEE;

  if(inputEpsilonFile_->IsOpen())
    inputEpsilonFile_->Close();
}

void FitEpsilonPlotForZ::loadEpsilonPlot(char *filename) {

  inputEpsilonFile_ = TFile::Open(filename);
  if(!inputEpsilonFile_) 
    throw cms::Exception("loadEpsilonPlot") << "Cannot open file " << string(filename) << "\n"; 

  char histoName[200];
  if( EEoEB_ == "Barrel" && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {

    for(int iR=inRangeFit_; iR <= finRangeFit_ && iR < nRegionsEB_; iR++) {
      sprintf(histoName, "WeightedRescaleFactorEB_channel_%d",iR);
      weightedRescaleFactorEB[iR] = (TH1F*)inputEpsilonFile_->Get(histoName);
      if(!weightedRescaleFactorEB[iR])
	throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << string(histoName) << "\n";
      else if(!(iR%1000))
	cout << "FIT_EPSILON_FORZ: Epsilon distribution for EB region " << iR << " loaded" << endl;
    }

  } else if( EEoEB_ == "Endcap" && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {

    for(int jR=inRangeFit_; jR <= finRangeFit_ && jR<EEDetId::kSizeForDenseIndexing; jR++) {
      sprintf(histoName, "WeightedRescaleFactorEE_channel_%d",jR);   
      weightedRescaleFactorEE[jR] = (TH1F*)inputEpsilonFile_->Get(histoName);
      if(!weightedRescaleFactorEE[jR])
	throw cms::Exception("loadEpsilonPlot") << "Cannot load histogram " << string(histoName) << "\n";
      else if(!(jR%1000))
	cout << "FIT_EPSILON: Epsilon distribution for EE region " << jR << " loaded" << endl;
    }
  }
}

void FitEpsilonPlotForZ::saveCoefficients()  {

  /// output file
  char fileName[200];
  sprintf(fileName,"%s/%s", outputDir_.c_str(), outfilename_.c_str());
  outfile_ = new TFile(fileName,"RECREATE");
  cout << "FIT_EPSILON_FOR_Z: Saving Calibration Coefficients in " << endl;
  cout << "FIT_EPSILON_FOR_Z: " << string(fileName) << " ... ";
  if(!outfile_) throw cms::Exception("WritingOutputFile") 
		  << "It was no possible to create output file " << string(fileName) << "\n";
  outfile_->cd();

  // 2D calib maps
  TH2F* hmap_EB = new TH2F("calibMap_EB","EB calib coefficients: #eta on x, #phi on y",
			   2*EBDetId::MAX_IETA+1,-EBDetId::MAX_IETA-0.5,EBDetId::MAX_IETA+0.5,
			   EBDetId::MAX_IPHI, EBDetId::MIN_IPHI-0.5, EBDetId::MAX_IPHI+0.5 );
  hmap_EB->GetXaxis()->SetTitle("i#eta");
  hmap_EB->GetYaxis()->SetTitle("i#phi");
  TH2F* hmap_EEp = new TH2F("calibMap_EEp","EE+ calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEp->GetXaxis()->SetTitle("ix");
  hmap_EEp->GetYaxis()->SetTitle("iy");
  TH2F* hmap_EEm = new TH2F("calibMap_EEm","EE- calib coefficients",100,0.5,100.5,100,0.5,100.5);
  hmap_EEm->GetXaxis()->SetTitle("ix");
  hmap_EEm->GetYaxis()->SetTitle("iy");
  TH1F* hint = new TH1F("hint","Bin1: inRangeFit_ Bin2: finRangeFit_ Bin3: Barrel(0)/Endcap(1)",3,0.,3.);
  hint->SetBinContent(1,inRangeFit_);
  hint->SetBinContent(2,finRangeFit_);
  if( EEoEB_ == "Barrel" ) hint->SetBinContent(3,0);
  else                     hint->SetBinContent(3,1);
  hint->Write();

  /// filling Barrel Map
  for(int j=0; j<nRegionsEB_; ++j) {
    std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
      EBDetId ebid(*iid);
      int ix = ebid.ieta()+EBDetId::MAX_IETA+1;
      float coeffValue = regionalCalibration_->getCalibMap()->coeff(*iid) > 0. ? regionalCalibration_->getCalibMap()->coeff(*iid) : 1.;
      hmap_EB->SetBinContent( ix, ebid.iphi(), coeffValue );
    } // loop over DetId in regions
  }
  hmap_EB->SetMinimum(0.9);
  hmap_EB->SetStats(false);
  hmap_EB->Write();

  /// filling Endcap Maps
  for(int jR=0; jR < regionalCalibration_->getCalibMap()->getNRegionsEE(); jR++) {
    std::vector<DetId> ids =  regionalCalibration_->allDetIdsInEERegion(jR);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) { 
      EEDetId eeid(*iid);
      float coeffValue =  regionalCalibration_->getCalibMap()->coeff(*iid) > 0. ?  regionalCalibration_->getCalibMap()->coeff(*iid) : 1.;
      if(eeid.positiveZ())
	hmap_EEp->Fill(eeid.ix(), eeid.iy(), coeffValue); 
      else 
	hmap_EEm->Fill(eeid.ix(), eeid.iy(), coeffValue);
    }
  }
  hmap_EEp->SetMinimum(0.9);
  hmap_EEp->SetStats(false);
  hmap_EEp->Write();  
  hmap_EEm->SetMinimum(0.9);
  hmap_EEm->SetStats(false);
  hmap_EEm->Write();
  
  // Output tree
  uint32_t   rawId;
  int        hashedIndex;
  int        ieta;
  int        iphi;
  int        iSM;
  int        iMod;
  int        iTT;
  int        iTTeta;
  int        iTTphi;
  int        iter = currentIteration_;
  float      regCoeff;
  float      Signal;//#
  float      Backgr; 
  float      Chisqu; 
  float      Ndof; 
  float      fit_mean;
  float      fit_mean_err;
  float      fit_sigma;
  float      fit_Snorm;
  float      fit_Bnorm; 
  int ix;
  int iy;
  int zside;
  int sc; 
  int isc;
  int ic;
  int iquadrant;

  TTree* treeEB = new TTree("calibEB","Tree of EB Inter-calibration constants");
  TTree* treeEE = new TTree("calibEE","Tree of EE Inter-calibration constants");

  /// barrel
  treeEB->Branch("rawId",&rawId,"rawId/i");
  treeEB->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEB->Branch("ieta",&ieta,"ieta/I");
  treeEB->Branch("iphi",&iphi,"iphi/I");
  treeEB->Branch("iSM",&iSM,"iSM/I");
  treeEB->Branch("iMod",&iMod,"iMod/I");
  treeEB->Branch("iTT",&iTT,"iTT/I");
  treeEB->Branch("iTTeta",&iTTeta,"iTTeta/I");
  treeEB->Branch("iTTphi",&iTTphi,"iTTphi/I");
  treeEB->Branch("iter",&iter,"iter/I");
  treeEB->Branch("coeff",&regCoeff,"coeff/F");
  treeEB->Branch("Signal",&Signal,"Signal/F");//#
  treeEB->Branch("Backgr",&Backgr,"Backgr/F");
  treeEB->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEB->Branch("Ndof",&Ndof,"Ndof/F");
  treeEB->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEB->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEB->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  treeEB->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  treeEB->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");
  
  /// endcap
  treeEE->Branch("ix",&ix,"ix/I");
  treeEE->Branch("iy",&iy,"iy/I");
  treeEE->Branch("zside",&zside,"zside/I");
  treeEE->Branch("sc",&sc,"sc/I");
  treeEE->Branch("isc",&isc,"isc/I");
  treeEE->Branch("ic",&ic,"ic/I");
  treeEE->Branch("iquadrant",&iquadrant,"iquadrant/I");
  treeEE->Branch("hashedIndex",&hashedIndex,"hashedIndex/I");
  treeEE->Branch("iter",&iter,"iter/I");
  treeEE->Branch("coeff",&regCoeff,"coeff/F");
  treeEE->Branch("Signal",&Signal,"Signal/F");//#
  treeEE->Branch("Backgr",&Backgr,"Backgr/F");
  treeEE->Branch("Chisqu",&Chisqu,"Chisqu/F");
  treeEE->Branch("Ndof",&Ndof,"Ndof/F");
  treeEE->Branch("fit_mean",&fit_mean,"fit_mean/F");
  treeEE->Branch("fit_mean_err",&fit_mean_err,"fit_mean_err/F");
  treeEE->Branch("fit_sigma",&fit_sigma,"fit_sigma/F");
  treeEE->Branch("fit_Snorm",&fit_Snorm,"fit_Snorm/F");
  treeEE->Branch("fit_Bnorm",&fit_Bnorm,"fit_Bnorm/F");
  
  
  for(int iR=0; iR < nRegionsEB_; ++iR)  {
    std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(iR);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
      EBDetId ebid(*iid);
      hashedIndex = ebid.hashedIndex();
      ieta = ebid.ieta();
      iphi = ebid.iphi();
      iSM = ebid.ism();
      iMod = ebid.im();
      iTT  = ebid.tower().hashedIndex();
      iTTeta = ebid.tower_ieta();
      iTTphi = ebid.tower_iphi();
      Signal = EBmap_Signal[ebid.hashedIndex()];//#
      Backgr = EBmap_Backgr[ebid.hashedIndex()];
      Chisqu = EBmap_Chisqu[ebid.hashedIndex()];
      Ndof = EBmap_ndof[ebid.hashedIndex()];
      fit_mean     = EBmap_mean[ebid.hashedIndex()];
      fit_mean_err = EBmap_mean_err[ebid.hashedIndex()];
      fit_sigma  = EBmap_sigma[ebid.hashedIndex()];
      fit_Snorm  = EBmap_Snorm[ebid.hashedIndex()];
      fit_Bnorm  = EBmap_Bnorm[ebid.hashedIndex()];
      regCoeff = regionalCalibration_->getCalibMap()->coeff(*iid);
      treeEB->Fill();
    } // loop over DetId in regions
  } // loop over regions

  for(int jR=0; jR < nRegionsEE_; jR++) {
    std::vector<DetId> ids = regionalCalibration_->allDetIdsInEERegion(jR);
    for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
      EEDetId eeid(*iid);
      ix = eeid.ix();
      iy = eeid.iy();
      zside = eeid.zside();
      sc = eeid.sc();
      isc = eeid.isc();
      ic = eeid.ic();
      iquadrant = eeid.iquadrant();
      hashedIndex = eeid.hashedIndex();
      regCoeff = regionalCalibration_->getCalibMap()->coeff(*iid);
      Signal = EEmap_Signal[eeid.hashedIndex()];//#
      Backgr = EEmap_Backgr[eeid.hashedIndex()];
      Chisqu = EEmap_Chisqu[eeid.hashedIndex()];            
      Ndof = EEmap_ndof[eeid.hashedIndex()];            
      fit_mean     = EEmap_mean[eeid.hashedIndex()];
      fit_mean_err = EEmap_mean_err[eeid.hashedIndex()];
      fit_sigma  = EEmap_sigma[eeid.hashedIndex()];
      fit_Snorm  = EEmap_Snorm[eeid.hashedIndex()];
      fit_Bnorm  = EEmap_Bnorm[eeid.hashedIndex()];
      treeEE->Fill();
    }
  }

  treeEB->Write();
  treeEE->Write();
  
  outfile_->Write();
  outfile_->Close();
  cout << "FIT_EPSILON_FOR_Z:  done with outputs" << endl;
}

// ------------ method called for each event  ------------
void FitEpsilonPlotForZ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  TF1 ffit("gausa","gaus(0)+[3]*x+[4]",-0.5,0.5);
  ffit.SetParameters(100,0,0.1);
  ffit.SetParNames("Constant","Mean_value","Sigma","a","b");
  ffit.SetParLimits(3,-500,500);
  ffit.SetParLimits(2,0.05,0.22);

  cout << "FIT_EPSILON_FOR_Z: About to fit epsilon distributions" << endl; 

  /// compute average weight, eps, and update calib constant
  if( (EEoEB_ == "Barrel") && (Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) ) {

    for(uint32_t j= (uint32_t)inRangeFit_; j <= (uint32_t)finRangeFit_ && j < (uint32_t)nRegionsEB_; ++j) { 
      cout<<"FIT_EPSILON_FOR_Z: Fitting EB Cristal--> "<<j<<endl;
      if(!(j%1000)) cout << "FIT_EPSILON_FOR_Z: fitting EB region " << j << endl;

      float mean = 0.;
      
      if(!useMassInsteadOfEpsilon_ && weightedRescaleFactorEB[j]->Integral(weightedRescaleFactorEB[j]->GetNbinsX()*(1./6.),weightedRescaleFactorEB[j]->GetNbinsX()*0.5) > 20) {
	double Max = 0.;
	double Min = -0.5, bin = 0.0125;
	Max = Min+(bin*(double)weightedRescaleFactorEB[j]->GetMaximumBin());
	double Bound1 = -0.15, Bound2 = 0.25;
	if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.1;}
	if ( Max+Bound2 > 0.34  ){ Bound2 = 0.15;}
	if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.075;}
	if ( Max+Bound2 > 0.34  ){ Bound2 = 0.1;}
	if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.03;}
	if ( Max+Bound2 > 0.34  ){ Bound2 = 0.05;}
	if ( fabs(Max+Bound1) > 0.24  ){ Bound1 = -0.009;}
	if ( Max+Bound2 > 0.34  ){ Bound2 = 0.01;}

	weightedRescaleFactorEB[j]->Fit(&ffit,"qB","", Max+Bound1,Max+Bound2);
	if(ffit.GetNDF() != 0) {
	  double chi2 = ( ffit.GetChisquare()/ffit.GetNDF() );
	  if ( chi2  > 11 ){
	    ffit.SetParLimits(2,0.05,0.15);
	    ffit.SetParameters(100,0,0.1);
	    weightedRescaleFactorEB[j]->Fit(&ffit,"qB","", Max+Bound1,Max+Bound2);
	    chi2 = (ffit.GetChisquare()/ffit.GetNDF());
	    if ( chi2  < 11 ){   cout<<"Saved 1 Level!!"<<endl;  }
	    else{
	      ffit.SetParameters(100,0,0.1);
	      ffit.SetParLimits(2,0.05,0.1);
	      weightedRescaleFactorEB[j]->Fit(&ffit,"qB","",  Max+Bound1,Max+Bound2);
	      chi2 = (ffit.GetChisquare()/ffit.GetNDF());
	      if ( chi2  < 11 ){ cout<<"Saved 2 Level!!"<<endl; }
	      else{ cout<<"DAMN: High Chi square..."<<endl; }
	    }
	  }
	}
	else cout<<"DAMN: NDF == 0"<<endl;
	mean = ffit.GetParameter(1);

      } else if (useMassInsteadOfEpsilon_) { 

	int iMin = weightedRescaleFactorEB[j]->GetXaxis()->FindBin(71);
	int iMax = weightedRescaleFactorEB[j]->GetXaxis()->FindBin(109);
	double integral = weightedRescaleFactorEB[j]->Integral(iMin, iMax);

	if(integral>60.) {        
	  ZFitResult fitres = FitMassPeakRooFit( weightedRescaleFactorEB[j], 70., 110., j, ZEB, 0); 
	  RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
	  mean = mean_fitresult->getVal();
	  float r2 = mean/MZ;
	  r2 = r2*r2;
	  float thechi2 = fitres.chi2;
	  float thendf  = fitres.dof;
	  float thechi2ndf = thechi2/thendf;
	  //if( thechi2<5 ) 
	  if( thechi2ndf<5 ) 
	    mean = 0.5 * ( r2 - 1. );
	  else  
	    mean = 0.;
	} else {
	  mean = 0.;
	}
      }

      std::vector<DetId> ids = regionalCalibration_->allDetIdsInEBRegion(j);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
	regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(1.+mean);
      } // loop over DetId in regions
    } // loop over regions
  } // if you have to fit barrel

  /// loop over EE crystals
  if( (EEoEB_ == "Endcap") && (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ){

    for(int jR = inRangeFit_; jR <=finRangeFit_ && jR < nRegionsEE_; jR++) {
      cout << "FIT_EPSILON_FOR_Z: Fitting EE Cristal--> " << jR << endl;
      if(!(jR%1000))
	cout << "FIT_EPSILON_FOR_Z: fitting EE region " << jR << endl;
      float mean = 0.;

      if(!useMassInsteadOfEpsilon_ && weightedRescaleFactorEE[jR]->Integral(weightedRescaleFactorEE[jR]->GetNbinsX()*(1./6.),weightedRescaleFactorEE[jR]->GetNbinsX()*0.5) > 20) {
	TF1 *ffit = new TF1("gausa","gaus(0)+[3]*x+[4]",-0.5,0.5);
	ffit->SetParameters(100,0,0.1);
	ffit->SetParNames("Constant","Mean_value","Sigma","a","b");
	ffit->SetParLimits(0,0.,weightedRescaleFactorEE[jR]->GetEntries()*1.1);
	ffit->SetParLimits(3,-500,500);
	ffit->SetParLimits(2,0.05,0.3);
	double Max = 0.;
	double Min = -0.5, bin = 0.0125;
	Max = Min+(bin*(double)weightedRescaleFactorEE[jR]->GetMaximumBin());
	double Bound1 = -0.35, Bound2 = 0.35;
	if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.3;}
	if ( Max+Bound2 > 0.48  ){ Bound2 = 0.3;}
	if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.25;}
	if ( Max+Bound2 > 0.48  ){ Bound2 = 0.2;}
	if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.2;}
	if ( Max+Bound2 > 0.48  ){ Bound2 = 0.15;}
	if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.15;}
	if ( Max+Bound2 > 0.48  ){ Bound2 = 0.1;}
	if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.1;}
	if ( fabs(Max+Bound1) > 0.38  ){ Bound1 = -0.05;}
	weightedRescaleFactorEE[jR]->Fit(ffit,"qB","", Max+Bound1,Max+Bound2);
	if(ffit->GetNDF() != 0) {
	  double chi2 = ( ffit->GetChisquare()/ffit->GetNDF() );
	  if(chi2 > 11  ) { cout<<"DAMN:(EE) High Chi square..."<<endl; }
	}
	else cout<<"DAMN: NDF == 0"<<endl;
	mean = ffit->GetParameter(1);

      } else if(useMassInsteadOfEpsilon_) {  

	int iMin = weightedRescaleFactorEE[jR]->GetXaxis()->FindBin(71);
	int iMax = weightedRescaleFactorEE[jR]->GetXaxis()->FindBin(109);
	double integral = weightedRescaleFactorEE[jR]->Integral(iMin, iMax);
	if(integral>60.) {
	  ZFitResult fitres = FitMassPeakRooFit( weightedRescaleFactorEE[jR], 70., 110., jR, ZEE, 0);  
	  RooRealVar* mean_fitresult = (RooRealVar*)(((fitres.res)->floatParsFinal()).find("mean"));
	  mean = mean_fitresult->getVal();
	  float r2 = mean/MZ;
	  r2 = r2*r2;
	  float thechi2 = fitres.chi2;
	  float thendf  = fitres.dof;
	  float thechi2ndf = thechi2/thendf;
	  //if( thechi2<5 ) 
	  if( thechi2ndf<5 ) 
	    mean = 0.5 * ( r2 - 1. );
	  else  
	    mean = 0.;
	} else {
	  mean = 0.;
	}
      }

      std::vector<DetId> ids = regionalCalibration_->allDetIdsInEERegion(jR);
      for(std::vector<DetId>::const_iterator iid = ids.begin(); iid != ids.end(); ++iid) {
	regionalCalibration_->getCalibMap()->coeff(*iid) *= (mean==0.) ? 1. : 1./(1.+mean);
      }
    } //for EE
  }  // if you have to fit Endcap
}

ZFitResult FitEpsilonPlotForZ::FitMassPeakRooFit(TH1F* h, double xlo, double xhi,  uint32_t HistoIndex, FitMode mode, int niter) {

  // variables
  RooRealVar x("x","ee invariant mass",xlo, xhi, "GeV/c^2");                  

  // signal
  RooRealVar mean ("mean", "mean", 91.,  80., 95.);
  RooRealVar sigma("sigma","sigma", 1., 0.01, 5.);     
  RooRealVar alpha("alpha","alpha",1.,0.,2.);
  RooRealVar CBn("CBn","CBn",1.5,0.,10.);
  RooCBShape signal("signal","signal", x, mean, sigma, alpha, CBn);

  // background
  RooRealVar expShape("expShape","expShape",-1.,-2.,0.);
  RooExponential bkg("background","background", x, expShape);

  // yields
  RooRealVar Nsig("Nsig","Z yield",1000.,0.,1.e7);
  Nsig.setVal( h->GetSum()*0.1);
  RooRealVar Nbkg("Nbkg","background yield",1.e3,0.,1.e8);
  Nbkg.setVal( h->GetSum()*0.8 );

  // full model
  RooAbsPdf* model=0;
  RooAddPdf model1("model","sig+bkg",RooArgList(signal,bkg),RooArgList(Nsig,Nbkg));
  model = &model1;

  // Fit
  RooDataHist dh("dh","ee invariant mass",RooArgList(x),h);
  RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
  RooMinuit m(nll);
  m.setVerbose(kFALSE);

  m.migrad();

  RooFitResult* res = m.save() ;

  // chi2 and S/B
  RooChi2Var chi2("chi2","chi2 var",*model,dh, true);
  int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();
  x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
  RooAbsReal* integralSig = signal.createIntegral(x,NormSet(x),Range("sobRange"));
  RooAbsReal* integralBkg = bkg.createIntegral(x,NormSet(x),Range("sobRange"));
  float normSig = integralSig->getVal();
  float normBkg = integralBkg->getVal();

  ZFitResult zres;     
  zres.res    = res;
  zres.S      = normSig*Nsig.getVal();
  zres.Serr   = normSig*Nsig.getError();
  zres.B      = normBkg*Nbkg.getVal();
  zres.Berr   = normBkg*Nbkg.getError();
  zres.SoB    = zres.S/zres.B;
  zres.SoBerr = zres.SoB*sqrt( pow(zres.Serr/zres.S,2) + pow(zres.Berr/zres.B,2) ) ;
  zres.dof    = ndof;

  RooPlot*  xframe = x.frame(h->GetNbinsX());
  xframe->SetTitle(h->GetTitle());
  dh.plotOn(xframe);
  model->plotOn(xframe,Components(bkg),LineStyle(kDashed), LineColor(kRed));
  model->plotOn(xframe);

  xframe->Draw();
  zres.probchi2 = TMath::Prob(xframe->chiSquare(), ndof);
  zres.chi2 = xframe->chiSquare();
  cout << "FIT_EPSILON_FORZ: Nsig: " << Nsig.getVal()
       << " nsig 3sig: " << normSig*Nsig.getVal()
       << " nbkg 3sig: " << normBkg*Nbkg.getVal()
       << " S/B: " << zres.SoB << " +/- " << zres.SoBerr
       << " chi2: " << xframe->chiSquare()
       << " DOF: " << zres.dof
       << " prob(chi2): " << zres.probchi2
       << endl;

  if(mode==ZEB){        
    EBmap_Signal[HistoIndex]=zres.S;
    EBmap_Backgr[HistoIndex]=zres.B;
    EBmap_Chisqu[HistoIndex]=xframe->chiSquare();
    EBmap_ndof[HistoIndex]=ndof;
    EBmap_mean[HistoIndex]=mean.getVal();
    EBmap_mean_err[HistoIndex]=mean.getError();
    EBmap_sigma[HistoIndex]=sigma.getVal();
    EBmap_Snorm[HistoIndex]=normSig;
    EBmap_Bnorm[HistoIndex]=normBkg;
  }
  if(mode==ZEE){        
    EEmap_Signal[HistoIndex]=zres.S;
    EEmap_Backgr[HistoIndex]=zres.B;
    EEmap_Chisqu[HistoIndex]=xframe->chiSquare();
    EEmap_ndof[HistoIndex]=ndof;
    EEmap_mean[HistoIndex]=mean.getVal();
    EEmap_mean_err[HistoIndex]=mean.getError();
    EEmap_sigma[HistoIndex]=sigma.getVal();
    EEmap_Snorm[HistoIndex]=normSig;
    EEmap_Bnorm[HistoIndex]=normBkg;
  }

  TLatex lat;
  char line[300];
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.SetTextColor(1);

  float xmin(0.55), yhi(0.80), ypass(0.05);
  if(mode==EtaEB) yhi=0.30;
  sprintf(line,"Yield: %.0f #pm %.0f", Nsig.getVal(), Nsig.getError() );
  lat.DrawLatex(xmin,yhi, line);
  sprintf(line,"m_{ee}: %.2f #pm %.2f", mean.getVal()*1000., mean.getError()*1000. );
  lat.DrawLatex(xmin,yhi-ypass, line);
  sprintf(line,"#sigma: %.2f #pm %.2f (%.2f%s)", sigma.getVal()*1000., sigma.getError()*1000., sigma.getVal()*100./mean.getVal(), "%" );
  lat.DrawLatex(xmin,yhi-2.*ypass, line);
  sprintf(line,"S/B(3#sigma): %.2f", zres.SoB );
  lat.DrawLatex(xmin,yhi-3.*ypass, line);
  sprintf(line,"#Chi^2: %.2f", xframe->chiSquare()/zres.dof );
  lat.DrawLatex(xmin,yhi-4.*ypass, line);
  sprintf(line,"Attempt: %d", niter );        
  lat.DrawLatex(xmin,yhi-5.*ypass, line);

  ZFitResult fitres = zres;               
  if(mode==ZEB && xframe->chiSquare()>5  ){
    if(niter==0) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, mode, 1);   
    if(niter==1) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, mode, 2);
    if(niter==2) fitres = FitMassPeakRooFit( h, xlo, xhi, HistoIndex, mode, 3);
  }

  if(StoreForTest_ && niter==0){
    std::stringstream ind;
    ind << (int) HistoIndex;
    TString nameHistofit = "Fit_n_" + ind.str();
    xframe->SetName(nameHistofit.Data());
    outfileTEST_->cd();
    xframe->Write();
  }
  
  return fitres;
}

void FitEpsilonPlotForZ::beginJob() {
  if(StoreForTest_){
    outfileTEST_ = new TFile("/tmp/Fit_Stored.root","RECREATE");
    if(!outfileTEST_) cout<<"WARNING: file with fit not created."<<endl;
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void FitEpsilonPlotForZ::endJob() {
  saveCoefficients();
  if(StoreForTest_){
    cout<<"Fit stored in /tmp/Fit_Stored.root"<<endl;
    outfileTEST_->Write();
    outfileTEST_->Close();
  }
}

void FitEpsilonPlotForZ::beginRun(edm::Run const&, edm::EventSetup const&) { }
void FitEpsilonPlotForZ::endRun(edm::Run const&, edm::EventSetup const&) { }
void FitEpsilonPlotForZ::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }
void FitEpsilonPlotForZ::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { } 
void FitEpsilonPlotForZ::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(FitEpsilonPlotForZ);
