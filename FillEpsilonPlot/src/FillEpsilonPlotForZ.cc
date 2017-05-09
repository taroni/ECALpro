// system include files
#include <memory>

// user include files
#include "TFile.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/CaloRecHit/interface/CaloID.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "CalibCode/FillEpsilonPlot/interface/FillEpsilonPlotForZ.h"
#include "CalibCode/FillEpsilonPlot/interface/JSON.h"

#define MZ 91.1876
#define MIN_RESCALE -0.5
#define MAX_RESCALE 0.5

#define MAX_PU_REWEIGHT 60

//#define DEBUG

FillEpsilonPlotForZ::FillEpsilonPlotForZ(const edm::ParameterSet& iConfig)
{
#ifdef DEBUG      
  cout << "DEBUG: Constructor" << endl;
#endif

  /// parameters from python: collections
  EBRecHitCollectionToken_  = consumes<EBRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EBRecHitCollectionTag"));
  EERecHitCollectionToken_  = consumes<EERecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EERecHitCollectionTag"));

  EBSuperClusterCollectionToken_  = consumes<SuperClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EBSuperClusterCollectionTag"));
  EESuperClusterCollectionToken_  = consumes<SuperClusterCollection>(iConfig.getUntrackedParameter<edm::InputTag>("EESuperClusterCollectionTag"));

  ElectronCollectionToken_  = consumes<GsfElectronCollection>(iConfig.getUntrackedParameter<edm::InputTag>("ElectronCollectionTag"));

  mcProducer_ = iConfig.getUntrackedParameter<std::string>("mcProducer");
  vertexToken_ = consumes<VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VertexTag"));
  PileUpToken_ = consumes<View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<InputTag> ("PileUpTag"));

  /// parameters from python: selection
  requireOppositeCharge_ = iConfig.getUntrackedParameter<bool>("requireOppositeCharge",false);
  minInvMassCut_         = iConfig.getUntrackedParameter<double>("minInvMassCut", 70.);
  maxInvMassCut_         = iConfig.getUntrackedParameter<double>("maxInvMassCut", 110.);
  electronSelection_     = iConfig.getUntrackedParameter<int> ("electronSelection",-1);
  elePtCut_              = iConfig.getUntrackedParameter<double>("elePtCut");
  eleEtaCut_             = iConfig.getUntrackedParameter<double>("eleEtaCut");
  maxDReleSc_            = iConfig.getUntrackedParameter<double>("maxDReleSc");
  massMethod_            = iConfig.getUntrackedParameter<string>("ZCalib_InvMass");

  /// parameters from python: others
  isMC_                    = iConfig.getUntrackedParameter<bool>("isMC",false);   
  isCRAB_                  = iConfig.getUntrackedParameter<bool>("isCRAB",false);
  Barrel_orEndcap_         = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");
  outputDir_               = iConfig.getUntrackedParameter<std::string>("OutputDir");
  outfilename_             = iConfig.getUntrackedParameter<std::string>("OutputFile");
  calibMapPath_            = iConfig.getUntrackedParameter<std::string>("calibMapPath");
  currentIteration_        = iConfig.getUntrackedParameter<int>("CurrentIteration");
  SystOrNot_               = iConfig.getUntrackedParameter<double>("SystOrNot",0);
  jsonFile_                = iConfig.getUntrackedParameter<std::string>("JSONfile","");
  useMassInsteadOfEpsilon_ = iConfig.getUntrackedParameter<bool>("useMassInsteadOfEpsilon",false);
  puWFileName_             = iConfig.getUntrackedParameter<std::string>("puWFileName");

  /// Json file
  if( jsonFile_!="" ) myjson=new JSON( edm::FileInPath( jsonFile_.c_str() ).fullPath().c_str() );

  /// setting calibration type: for the moment per crystal only implemented
  calibTypeNumber_ = xtal;    
  regionalCalibration_ = &xtalCalib; 
#ifdef DEBUG
  cout << "crosscheck: selected type: " << regionalCalibration_->printType() << endl;
#endif

  /// retrieving calibration coefficients of the previous iteration
  if(currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
  else if(currentIteration_ > 0 || calibMapPath_.find("iter_-1")==std::string::npos) {
    char fileName[200];
    cout << "FillEpsilonPlotForZ:: loading calibraion map at " << calibMapPath_ << endl;
    if( isCRAB_ ) sprintf(fileName,"%s",  edm::FileInPath( calibMapPath_.c_str() ).fullPath().c_str() );
    else          sprintf(fileName,"%s", calibMapPath_.c_str());
    regionalCalibration_->getCalibMap()->loadCalibMapFromFile(fileName);
  }

  // # of regions
  nRegionsEB_ = regionalCalibration_->getCalibMap()->getNRegionsEB() ;
  nRegionsEE_ = regionalCalibration_->getCalibMap()->getNRegionsEE() ;
  if (nRegionsEB_>N_XTAL_BARREL) throw cms::Exception("WritingOutputFile") << "FillEpsilonPlotForZ => wring number of EB regions \n";
  if (nRegionsEE_>N_XTAL_ENDCAP) throw cms::Exception("WritingOutputFile") << "FillEpsilonPlotForZ => wring number of EE regions \n";

  // Histograms
  bookHistograms();

  /// Output file
  char fileName[200];
  sprintf(fileName,"%s%s", outputDir_.c_str(), outfilename_.c_str());
  outfile_ = new TFile(fileName,"RECREATE");
  if(!outfile_) throw cms::Exception("WritingOutputFile") << "It was no possible to create output file " << string(fileName) << "\n";

  /// Output tree
  myTree = new TTree("myTree","myTree");
  myTree->Branch("iteration",&currentIteration_,"iteration/I");
  myTree->Branch("isEBEB",&isEBEB,"isEBEB/I");
  myTree->Branch("isEEEE",&isEEEE,"isEEEE/I");
  myTree->Branch("isHR9HR9",&isHR9HR9,"isHR9HR9/I");
  myTree->Branch("nvtx",&nvtx,"nvtx/I");
  myTree->Branch("zMass",&mass4tree,"mass/F");
  myTree->Branch("mcweight",&thisEventW,"mcweight/F");
  myTree->Branch("puweight",&pu_weight,"puweight/F");
  myTree->Branch("ptele1",&ptele1,"ptele1/F");
  myTree->Branch("ptele2",&ptele2,"ptele2/F");
  myTree->Branch("etaele1",&etaele1,"etaele1/F");
  myTree->Branch("etaele2",&etaele2,"etaele2/F");

#ifdef DEBUG      
  cout << "DEBUG: Constructor done" << endl;
#endif
}

FillEpsilonPlotForZ::~FillEpsilonPlotForZ() {

#ifdef DEBUG      
  cout << "DEBUG: Destructor" << endl;
#endif

  outfile_->Write();
  outfile_->Close();

  if( Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) deleteEpsilonPlot(weightedRescaleFactorEB, nRegionsEB_);
  if( Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) deleteEpsilonPlot(weightedRescaleFactorEE, nRegionsEE_);

  delete EventFlow;
  delete allEpsilon_EB;
  delete allEpsilon_EBnw;
  delete allEpsilon_EE;
  delete allEpsilon_EEnw;

  delete myTree;

  //JSON
  if( jsonFile_!="" ) delete myjson;

#ifdef DEBUG      
  cout << "DEBUG: Destructor" << endl;
#endif
}

// ------------ method called for each event  ------------
void FillEpsilonPlotForZ::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
#ifdef DEBUG      
  cout << "DEBUG: Analyze" << endl;
#endif

  
  // PU weight (for MC only)
  pu_weight = 1.; 
  float pu_n = -1.;
  if (isMC_) {  
    Handle<View< PileupSummaryInfo> > PileupInfos;
    iEvent.getByToken(PileUpToken_,PileupInfos);
    pu_n = 0.;
    for( unsigned int PVI = 0; PVI < PileupInfos->size(); ++PVI ) {
      Int_t pu_bunchcrossing = PileupInfos->ptrAt( PVI )->getBunchCrossing();
      if( pu_bunchcrossing == 0 ) 
	pu_n = PileupInfos->ptrAt( PVI )->getTrueNumInteractions();         
    }
    pu_weight = GetPUWeight(pu_n);         
  }

  // dataset weight (for Mc only)
  float perEveW = 1.;
  if (isMC_ && !mcProducer_.empty()) {
    Handle<GenEventInfoProduct> genInfo;
    try {
      iEvent.getByLabel("generator",genInfo);
    } catch (std::exception& ex) {
      std::cerr << "Error! can't get the product genInfo" << std::endl;
    }
    const auto & eveWeights = genInfo->weights();
    if(!eveWeights.empty()) perEveW = eveWeights[0];
  }

  // Total MC weight
  thisEventW = pu_weight*perEveW;

  // All events
  EventFlow->Fill(0.,thisEventW); 

  // Json file
  if ( jsonFile_!="" && !myjson->isGoodLS(iEvent.id().run(),iEvent.id().luminosityBlock()) ) return;
  EventFlow->Fill(1.,thisEventW); 

  // For Syst error SystOrNot_=1 or 2, for normal calib is 0
  if(SystOrNot_==1. && int(iEvent.id().event())%2!=0 )      return;
  else if(SystOrNot_==2. && int(iEvent.id().event())%2==0 ) return;


  // ----------------------------------------------------------
  // Loading collections

  // vertices
  iEvent.getByToken( vertexToken_, primaryVertices );
  nvtx = primaryVertices->size();

  // ECAL Rechits
  iEvent.getByToken ( EBRecHitCollectionToken_, hits);
  iEvent.getByToken ( EERecHitCollectionToken_, ehits);
  if (hits->size()==0 && ehits->size()==0) {
    std::cout << "hits->size() == 0" << std::endl;
    return;
  }

  // SuperClusters
  iEvent.getByToken ( EBSuperClusterCollectionToken_, ebScCollection );
  iEvent.getByToken ( EESuperClusterCollectionToken_, eeScCollection );
  if (ebScCollection->size()==0 && eeScCollection->size()==0) {
    std::cout << "superclusters->size() == 0" << std::endl;
    return;
  }

#ifdef DEBUG      
  cout << "DEBUG: ";
  cout << "taken: "
       << ebScCollection->size() << " sc in EB, "
       << eeScCollection->size() << " sc in EE, "
       << hits->size()           << " rechits in EB, "
       << ehits->size()          << " rechits in EE " << endl;
#endif

#ifdef DEBUG_COLL
  cout << "DEBUG: ";
  std::cout<<"EB: ebScCollection->size() = " << ebScCollection->size() << std::endl;
  for(reco::SuperClusterCollection::const_iterator scIt = ebScCollection->begin(); scIt != ebScCollection->end(); scIt++) {
    std::cout << scIt->energy() << " " << scIt->eta() << " " << scIt->phi() << std::endl;
  }
  std::cout<<"EE: eeScCollection->size() = " << eeScCollection->size() << std::endl;
  for(reco::SuperClusterCollection::const_iterator scIt = eeScCollection->begin(); scIt != eeScCollection->end(); scIt++) {
    std::cout << scIt->energy() << " " << scIt->eta() << " " << scIt->phi() << std::endl;
  }
  std::cout << std::endl;
#endif

  if(  ( ebScCollection->size()+eeScCollection->size() ) < 2) return;
  EventFlow->Fill(2.,thisEventW); 

  // Electrons
  iEvent.getByToken ( ElectronCollectionToken_, electronCollection );
  
#ifdef DEBUG
  cout << "DEBUG: ";
  std::cout<<"electronCollection->size() = " << electronCollection->size() << std::endl;
  for(reco::GsfElectronCollection::const_iterator eleIt = electronCollection->begin(); eleIt != electronCollection->end(); eleIt++) {
    std::cout << eleIt->energy() << " " << eleIt->eta() << " " << eleIt->phi() << std::endl;
  }
  std::cout << std::endl;
#endif
  
  if(electronCollection->size() < 2) return;
  EventFlow->Fill(3.,thisEventW);


  // ----------------------------------------------------------
  // Building electrons to make Zs, using recalibrated SC
  std::vector<calib::CalibElectronForZ> calibElectrons;         

  // Geometrical match between (recalibrated) superclusters and (standard) electrons 
  for(reco::GsfElectronCollection::const_iterator eleIt = electronCollection->begin(); eleIt != electronCollection->end(); eleIt++) {

    if (!(eleIt->parentSuperCluster())) continue;

    // only electrons in the acceptance and with pT above threshold 
    float etaFromSc   = eleIt->parentSuperCluster()->eta();
    float thetaFromSc = eleIt->parentSuperCluster()->position().theta();
    float ptFromSc    = (eleIt->parentSuperCluster()->energy())*fabs(sin(thetaFromSc));
    if (ptFromSc<elePtCut_)         continue;               
    if (fabs(etaFromSc)>eleEtaCut_) continue;

    // ele-SC matching -    // chiara: gap EB/EE excluded                                                          
    float DeltaRMineleSCbarrel(maxDReleSc_);
    float DeltaRMineleSCendcap(maxDReleSc_);

    const reco::SuperCluster* sc=0;
    if(eleIt->isEB()) {
      for(reco::SuperClusterCollection::const_iterator scEbIt = ebScCollection->begin(); scEbIt != ebScCollection->end(); scEbIt++) {
        double DeltaReleSC = deltaR(eleIt->parentSuperCluster()->eta(),eleIt->parentSuperCluster()->phi(),scEbIt->eta(),scEbIt->phi());
        if(DeltaReleSC<DeltaRMineleSCbarrel) {
          DeltaRMineleSCbarrel = DeltaReleSC;
          sc=&(*scEbIt);
        }
      }
    }
    else if(eleIt->isEE()) {
      for(reco::SuperClusterCollection::const_iterator scEeIt = eeScCollection->begin(); scEeIt != eeScCollection->end(); scEeIt++) {
        double DeltaReleSC = deltaR(eleIt->parentSuperCluster()->eta(),eleIt->parentSuperCluster()->phi(),scEeIt->eta(),scEeIt->phi());
        if(DeltaReleSC<DeltaRMineleSCendcap) {
          DeltaRMineleSCendcap = DeltaReleSC;
          sc = &(*scEeIt);
        }
      }
    }
    if(!sc) continue;
    
    calibElectrons.push_back(calib::CalibElectronForZ(&(*eleIt),sc,&(*hits),&(*ehits)));

#ifdef DEBUG
    cout << "DEBUG: ";
    std::cout << "calibElectrons.back().getParentSuperCluster()->energy() = "
              << calibElectrons.back().getParentSuperCluster()->energy()
              << ", calibElectrons.back().getRecoElectron()->energy() = "
              << calibElectrons.back().getRecoElectron()->energy() << std::endl;
#endif
  }
  
#ifdef DEBUG
  cout << "DEBUG: ";
  std::cout << "calibElectrons.size() = " << calibElectrons.size() << endl;
#endif

#ifdef DEBUG
  cout << "DEBUG: ";
  std::cout << "Filled histos" << std::endl;
#endif

  if (calibElectrons.size() < 2) return;
  EventFlow->Fill(4.,thisEventW);


  // ----------------------------------------------------------
  // Combinatorics for Z mass
  std::vector<std::pair<calib::CalibElectronForZ*,calib::CalibElectronForZ*> > zeeCandidates;
  int  myBestZ = -1;
  float tmpMass = -1.;
  double DeltaMinvMin(5000.);

  for(unsigned int e_it = 0 ; e_it != calibElectrons.size() - 1 ; e_it++) {
    for(unsigned int p_it = e_it + 1 ; p_it != calibElectrons.size() ; p_it++) {

#ifdef DEBUG
      cout << "DEBUG: ";
      std::cout << e_it << " " << calibElectrons[e_it].getRecoElectron()->charge() << " "
                << p_it << " " << calibElectrons[p_it].getRecoElectron()->charge() << std::endl;
#endif

      // opposite charge   
      if ( requireOppositeCharge_ && (calibElectrons[e_it].getRecoElectron()->charge() * calibElectrons[p_it].getRecoElectron()->charge() != -1) ) continue;

      // when selecting the same SC for the two electrons I drop the event
      if (calibElectrons[e_it].getParentSuperCluster() == calibElectrons[p_it].getParentSuperCluster()) continue;

      // invariant mass
      if (massMethod_ == "SCTRMass" ) {
	tmpMass = invMassCalc(calibElectrons[p_it].getParentSuperCluster()->energy(), calibElectrons[p_it].getRecoElectron()->eta(), calibElectrons[p_it].getRecoElectron()->phi(), calibElectrons[e_it].getParentSuperCluster()->energy(), calibElectrons[e_it].getRecoElectron()->eta(), calibElectrons[e_it].getRecoElectron()->phi());
      } else if (massMethod_ == "SCMass" ) {
	tmpMass = invMassCalc(calibElectrons[p_it].getParentSuperCluster()->energy(), calibElectrons[p_it].getParentSuperCluster()->position().eta(), calibElectrons[p_it].getParentSuperCluster()->position().phi(), calibElectrons[e_it].getParentSuperCluster()->energy(), calibElectrons[e_it].getParentSuperCluster()->position().eta(), calibElectrons[e_it].getParentSuperCluster()->position().phi());
      }  
      
      if (tmpMass<0) continue;

#ifdef DEBUG
      cout << "DEBUG: ";
      std::cout << "####################### mass = " << tmpMass << std::endl;
#endif
 
      zeeCandidates.push_back(std::pair<calib::CalibElectronForZ*,calib::CalibElectronForZ*>(&(calibElectrons[e_it]),&(calibElectrons[p_it])));

      double DeltaMinv = fabs(tmpMass - MZ);
      if( DeltaMinv < DeltaMinvMin) {
        DeltaMinvMin = DeltaMinv;
        myBestZ=zeeCandidates.size()-1;
      }
    }
  }

  if(zeeCandidates.size()==0 || myBestZ==-1 ) return;
  EventFlow->Fill(5.,thisEventW);

#ifdef DEBUG
  cout << "DEBUG: ";
  std::cout << "Found ZCandidates: " << myBestZ << std::endl;
#endif

  // ----------------------------------------------------------
  // Mass window
  bool invMassBool = ( (tmpMass > minInvMassCut_) && (tmpMass < maxInvMassCut_) );
  if (!invMassBool) return;
  EventFlow->Fill(6.,thisEventW);

#ifdef DEBUG
  cout << "DEBUG: ";
  std::cout << "Found ZCandidates within mass window: " << myBestZ << std::endl;
#endif

  // ----------------------------------------------------------
  // Classification:
  // enum Classification { UNKNOWN=-1, GOLDEN=0, BIGBREM=1, BADTRACK=2, SHOWERING=3, GAP=4 } 
  // in DataFormats/EgammaCandidates/interface/GsfElecron.h

  float eta1 = zeeCandidates[myBestZ].first->getParentSuperCluster()->eta();
  float eta2 = zeeCandidates[myBestZ].second->getParentSuperCluster()->eta();
  
  float r91  = zeeCandidates[myBestZ].first->getRecoElectron()->r9(); 
  float r92  = zeeCandidates[myBestZ].second->getRecoElectron()->r9(); 

  bool selectionBool = false;  
  // -1 = all electrons
  if(electronSelection_==-1)     selectionBool=( myBestZ != -1 );    
  // 0 = all electrons but not gap
  else if(electronSelection_==0) selectionBool=( myBestZ != -1 && 
						 zeeCandidates[myBestZ].first->getRecoElectron()->classification()!= 4 && 
						 zeeCandidates[myBestZ].second->getRecoElectron()->classification()!= 4 ); 
  if (!selectionBool) return;
  EventFlow->Fill(7.,thisEventW);

  // reject crazy electrons
  if ( zeeCandidates[myBestZ].first->getParentSuperCluster()->position().eta() < -10.  ) return; 
  if ( zeeCandidates[myBestZ].second->getParentSuperCluster()->position().eta() < -10. ) return; 
  if ( zeeCandidates[myBestZ].first->getParentSuperCluster()->position().phi() < -10.  ) return; 
  if ( zeeCandidates[myBestZ].second->getParentSuperCluster()->position().phi() < -10. ) return; 
  EventFlow->Fill(8.,thisEventW);

  // TLorentzVectors
  //math::PtEtaPhiMLorentzVector e1P4( (zeeCandidates[myBestZ].first->getParentSuperCluster()->energy()/zeeCandidates[myBestZ].first->getRecoElectron()->eta()), zeeCandidates[myBestZ].first->getRecoElectron()->eta(), zeeCandidates[myBestZ].first->getRecoElectron()->phi(), 0. );
  //math::PtEtaPhiMLorentzVector e2P4( (zeeCandidates[myBestZ].second->getParentSuperCluster()->energy()/zeeCandidates[myBestZ].second->getRecoElectron()->eta()), zeeCandidates[myBestZ].second->getRecoElectron()->eta(), zeeCandidates[myBestZ].second->getRecoElectron()->phi(), 0. );
  //math::PtEtaPhiMLorentzVector zP4 = e1P4 + e2P4;
  

  // -----------------------------------------------------------
  // Real calibration step
  if (massMethod_ == "SCTRMass" ) {
    mass4tree = invMassCalc(zeeCandidates[myBestZ].first->getParentSuperCluster()->energy(), zeeCandidates[myBestZ].first->getRecoElectron()->eta(), zeeCandidates[myBestZ].first->getRecoElectron()->phi(), zeeCandidates[myBestZ].second->getParentSuperCluster()->energy(), zeeCandidates[myBestZ].second->getRecoElectron()->eta(), zeeCandidates[myBestZ].second->getRecoElectron()->phi());
  } else if (massMethod_ == "SCMass" ) {
    mass4tree = invMassCalc(zeeCandidates[myBestZ].first->getParentSuperCluster()->energy(), zeeCandidates[myBestZ].first->getParentSuperCluster()->position().eta(), zeeCandidates[myBestZ].first->getParentSuperCluster()->position().phi(), zeeCandidates[myBestZ].second->getParentSuperCluster()->energy(), zeeCandidates[myBestZ].second->getParentSuperCluster()->position().eta(), zeeCandidates[myBestZ].second->getParentSuperCluster()->position().phi());
  }  

  if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) {
    if (fabs(eta1)<1.5) getWeight(mass4tree, zeeCandidates[myBestZ].first,  MZ, thisEventW, EcalBarrel, &(*hits));     
    if (fabs(eta2)<1.5) getWeight(mass4tree, zeeCandidates[myBestZ].second, MZ, thisEventW, EcalBarrel, &(*hits));     
  }
  if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) {
    if (fabs(eta1)>1.5) getWeight(mass4tree, zeeCandidates[myBestZ].first,  MZ, thisEventW, EcalEndcap, &(*ehits));     
    if (fabs(eta2)>1.5) getWeight(mass4tree, zeeCandidates[myBestZ].second, MZ, thisEventW, EcalEndcap, &(*ehits));     
  }

  // -----------------------------------------------------------
  // Tree filling 
  isEBEB=0;
  isEEEE=0;
  if (fabs(eta1)<1.5 && fabs(eta2)<1.5 ) isEBEB=1;
  if (fabs(eta1)>1.5 && fabs(eta2)>1.5 ) isEEEE=1;

  isHR9HR9=0;
  if (r91>0.94 && r92>0.94) isHR9HR9 = 1;

  ptele1 = zeeCandidates[myBestZ].first->getParentSuperCluster()->energy()*sin(zeeCandidates[myBestZ].first->getRecoElectron()->theta());
  ptele2 = zeeCandidates[myBestZ].second->getParentSuperCluster()->energy()*sin(zeeCandidates[myBestZ].second->getRecoElectron()->theta());
  etaele1 = eta1;
  etaele2 = eta2;

  myTree->Fill();
}

// filled with
// mreco(for this event) 
// OR
// [( mreco(for this event) / mZ(from PDG, with corrections) )^2 -1] / 2
// 
// weighted with
// weight2 = rawEne of SC / sum of raw energies of the SC hits in each module

void FillEpsilonPlotForZ::getWeight(float recomass, calib::CalibElectronForZ* ele, float evweight, float mcWeight, int subDetId, const EcalRecHitCollection* theHits) {
  
  RegionWeightVector wEle = regionalCalibration_->getWeightsZ( ele->getParentSuperCluster(), subDetId, theHits );     

  float tofill = recomass;

  if (!useMassInsteadOfEpsilon_) tofill = (TMath::Power((recomass / evweight), 2.) - 1) / 2.;
  if (!useMassInsteadOfEpsilon_ && (tofill< MIN_RESCALE || tofill>MAX_RESCALE)) {
    std::cout << "[FillEpsilonPlotForZ]::[getWeight]::rescale out " << tofill << std::endl;

  } else {

    if(subDetId!=EcalBarrel) allEpsilon_EEnw->Fill( tofill );
    if(subDetId==EcalBarrel) allEpsilon_EBnw->Fill( tofill );

    for(RegionWeightVector::const_iterator it = wEle.begin(); it != wEle.end(); ++it) {

      float corrToRawWeight = ele->getParentSuperCluster()->energy()/ele->getParentSuperCluster()->rawEnergy();

      const uint32_t& iR   = (*it).iRegion;
      const float& weight2 = (*it).value * corrToRawWeight;

      if (weight2>=0. && weight2<=1.) {
	
	if (subDetId==EcalBarrel) {
	  weightedRescaleFactorEB[iR]->Fill(tofill,weight2*mcWeight);
	  allEpsilon_EB->Fill(tofill,weight2*mcWeight);
	} else if (subDetId==EcalEndcap) {
	  allEpsilon_EE->Fill(tofill,weight2*mcWeight);  
	  weightedRescaleFactorEE[iR]->Fill(tofill,weight2*mcWeight);
	}
      }
    }
  }

}

// ------------ method called once each job just before starting event loop  ------------
void FillEpsilonPlotForZ::beginJob()
{
#ifdef DEBUG
  cout << "DEBUG: ";
  cout << "[DEBUG] beginJob" << endl;
#endif

  // loading pileup weights
  if (isMC_) SetPuWeights(puWFileName_);    
}

void FillEpsilonPlotForZ::endJob(){

  outfile_->cd();
  myTree->Write();

  EventFlow->Write();
  allEpsilon_EB->Write();
  allEpsilon_EBnw->Write();
  allEpsilon_EE->Write();
  allEpsilon_EEnw->Write();

  if( Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) 
    writeEpsilonPlot(weightedRescaleFactorEB,nRegionsEB_);

  if( Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) 
    writeEpsilonPlot(weightedRescaleFactorEE,nRegionsEE_);
}

void FillEpsilonPlotForZ::SetPuWeights(std::string puWeightFile) {

  if (puWeightFile == "") {
    std::cout << "you need a weights file to use this function" << std::endl;
    return;
  }
  std::cout << "PU REWEIGHTING:: Using file " << puWeightFile << std::endl;

  TFile *f_pu  = new TFile(puWeightFile.c_str(),"READ");
  f_pu->cd();
  
  TH1D *puweights = 0;
  TH1D *gen_pu = 0;
  gen_pu    = (TH1D*) f_pu->Get("generated_pu");
  puweights = (TH1D*) f_pu->Get("weights");
  
  if (!puweights || !gen_pu) {
    std::cout << "weights histograms  not found in file " << puWeightFile << std::endl;
    return;
  }
  TH1D* weightedPU= (TH1D*)gen_pu->Clone("weightedPU");
  weightedPU->Multiply(puweights);
  
  // Rescaling weights in order to preserve same integral of events                               
  TH1D* weights = (TH1D*)puweights->Clone("rescaledWeights");
  weights->Scale( gen_pu->Integral(1,MAX_PU_REWEIGHT) / weightedPU->Integral(1,MAX_PU_REWEIGHT) );
  
  float sumPuWeights=0.;
  for (int i = 0; i<MAX_PU_REWEIGHT; i++) {
    float weight=1.;
    weight=weights->GetBinContent(i+1);
    sumPuWeights+=weight;
    puweights_.push_back(weight);
  }
}

float FillEpsilonPlotForZ::GetPUWeight(float pun) {
  
  float weight=1;
  if (isMC_ && pun<MAX_PU_REWEIGHT && puweights_.size()>0) 
    weight = puweights_[pun];
  return weight;
}

// histo booking
void FillEpsilonPlotForZ::bookHistograms() {

  int nbins = 100;
  float lowH  = MIN_RESCALE;
  float highH = MAX_RESCALE;
  if (useMassInsteadOfEpsilon_) {
    nbins = maxInvMassCut_-minInvMassCut_;
    lowH  = minInvMassCut_;
    highH = maxInvMassCut_;
  }

  if( Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) 
    weightedRescaleFactorEB   = initializeEpsilonHistograms("WeightedRescaleFactorEB_channel_",  "WeightedRescaleFactorEB_channel_",   nRegionsEB_,0);

  if( (Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) ) 
    weightedRescaleFactorEE   = initializeEpsilonHistograms("WeightedRescaleFactorEE_channel_",  "WeightedRescaleFactorEE_channel_",   nRegionsEE_,0);

  EventFlow = new TH1F("EventFlow", "EventFlow", 10, -0.5, 9.5 );
  EventFlow->GetXaxis()->SetBinLabel(1,"All Events"); 
  EventFlow->GetXaxis()->SetBinLabel(2,"JSON"); 
  EventFlow->GetXaxis()->SetBinLabel(3,"Existing SCs"); 
  EventFlow->GetXaxis()->SetBinLabel(4,"Existing Eles"); 
  EventFlow->GetXaxis()->SetBinLabel(5,"Match Ele-SCs"); 
  EventFlow->GetXaxis()->SetBinLabel(6,"Existing Zs"); 
  EventFlow->GetXaxis()->SetBinLabel(7,"m(ee) cut"); 
  EventFlow->GetXaxis()->SetBinLabel(8,"Ele selection"); 
  EventFlow->GetXaxis()->SetBinLabel(9,"Meaningfull Eles"); 

  allEpsilon_EB   = new TH1F("allEpsilon_EB",   "allEpsilon_EB",   nbins, lowH, highH);
  allEpsilon_EBnw = new TH1F("allEpsilon_EBnw", "allEpsilon_EBnw", nbins, lowH, highH);
  allEpsilon_EE   = new TH1F("allEpsilon_EE",   "allEpsilon_EE",   nbins, lowH, highH);
  allEpsilon_EEnw = new TH1F("allEpsilon_EEnw", "allEpsilon_EEnw", nbins, lowH, highH);
}

void FillEpsilonPlotForZ::deleteEpsilonPlot(TH1F **h, int size) {
  for(int jR=0; jR<size; jR++)
    delete h[jR];

  delete h;
}

void FillEpsilonPlotForZ::writeEpsilonPlot(TH1F **h, int size)
{
  for(int jR=0; jR<size; jR++) h[jR]->Write();
}

TH1F** FillEpsilonPlotForZ::initializeEpsilonHistograms(const char *name, const char *title, int size, int isweight ) {

  TH1F **h = new TH1F*[size];
  char name_c[200];
  char title_c[200];
  
  int nbins = 100;
  float lowH  = MIN_RESCALE;
  float highH = MAX_RESCALE;
  if (useMassInsteadOfEpsilon_) {
    nbins = maxInvMassCut_-minInvMassCut_;
    lowH  = minInvMassCut_;
    highH = maxInvMassCut_;
  }
  
  for(int jR=0; jR<size; jR++){
    sprintf(name_c,  "%s%d", name, jR);
    sprintf(title_c, "%s%d", title, jR);
    h[jR] = new TH1F(name_c, title_c, nbins, lowH, highH);
    if (!isweight) {
      if(useMassInsteadOfEpsilon_) h[jR]->GetXaxis()->SetTitle("Z mass");
      else h[jR]->GetXaxis()->SetTitle("Rescale factor");    
    } else {
      h[jR]->GetXaxis()->SetTitle("Weight");
    }
  }

  return h;
}

void FillEpsilonPlotForZ::beginRun(edm::Run const&, edm::EventSetup const& iSetup) { }

void FillEpsilonPlotForZ::endRun(edm::Run const&, edm::EventSetup const&) { }

void FillEpsilonPlotForZ::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

void FillEpsilonPlotForZ::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) { }

void FillEpsilonPlotForZ::fillDescriptions(edm::ConfigurationDescriptions& descriptions) { }     

//define this as a plug-in
DEFINE_FWK_MODULE(FillEpsilonPlotForZ);
