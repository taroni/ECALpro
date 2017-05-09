#include "CalibCode/FillEpsilonPlot/plugins/CalibRechitsForZProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include <iostream>

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

using namespace edm;
using namespace std;

CalibRechitsForZProducer::CalibRechitsForZProducer(const edm::ParameterSet& iConfig)
{
  /// Parameters relevant for the calibration
  currentIteration_ = iConfig.getUntrackedParameter<int>("CurrentIteration");
  calibMapPath_     = iConfig.getUntrackedParameter<std::string>("calibMapPath");
  isCRAB_           = iConfig.getUntrackedParameter<bool>("isCRAB",false);
  Barrel_orEndcap_  = iConfig.getUntrackedParameter<std::string>("Barrel_orEndcap");

  /// Input collections
  EBRecHitCollectionToken_  = consumes<EBRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollectionTag"));
  EERecHitCollectionToken_  = consumes<EERecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHitCollectionTag"));

  /// Output collections
  EBNewRecHitCollection_ = iConfig.getParameter<std::string>("EBNewRecHitCollection");
  EENewRecHitCollection_ = iConfig.getParameter<std::string>("EENewRecHitCollection");

  /// Setting calibration type - for the moment, crystal granularity only implemented
  regionalCalibration_ = &xtalCalib; 
  cout << "crosscheck: selected type: " << regionalCalibration_->printType() << endl;

  /// Retrieving calibration coefficients of the previous iteration
  if(currentIteration_ < 0) throw cms::Exception("IterationNumber") << "Invalid negative iteration number\n";
  else if(currentIteration_ > 0 || calibMapPath_.find("iter_-1")==std::string::npos) {
    char fileName[200];
    cout << "CalibRechitsForZ:: loading calibraion map at " << calibMapPath_ << endl;
    if( isCRAB_ ) sprintf(fileName,"%s",  edm::FileInPath( calibMapPath_.c_str() ).fullPath().c_str() );
    else          sprintf(fileName,"%s", calibMapPath_.c_str());
    regionalCalibration_->getCalibMap()->loadCalibMapFromFile(fileName);
  }

  /// Outputs setting
  produces< EBRecHitCollection >(EBNewRecHitCollection_);
  produces< EERecHitCollection >(EENewRecHitCollection_);
}

CalibRechitsForZProducer::~CalibRechitsForZProducer() { }

void CalibRechitsForZProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // reading rechits before previous round calibration
  iEvent.getByToken ( EBRecHitCollectionToken_, ebHandle);
  iEvent.getByToken ( EERecHitCollectionToken_, eeHandle);
  
  // collections to be put in the event
  std::auto_ptr< EBRecHitCollection > EBNewRecHits( new EBRecHitCollection );
  std::auto_ptr< EERecHitCollection > EENewRecHits( new EERecHitCollection );

  // loop over barrel rechits if wanted
  if(Barrel_orEndcap_=="ONLY_BARREL" || Barrel_orEndcap_=="ALL_PLEASE" ) {
   
    for(EBRecHitCollection::const_iterator itb= ebHandle->begin(); itb != ebHandle->end(); ++itb ) { 

      EBDetId tmp_id(itb->id());
      float rh_calib = regionalCalibration_->getCalibMap()->coeff(tmp_id); 

      // calibrated energy
      float tmp_ene = itb->energy();
      float calib_ene = tmp_ene * rh_calib;

      // new rechit
      EcalRecHit aHit( *itb );
      aHit.setEnergy(calib_ene);
      EBNewRecHits->push_back( aHit );
    }
  }

  // loop over endcap rechits if wanted
  if(Barrel_orEndcap_=="ONLY_ENDCAP" || Barrel_orEndcap_=="ALL_PLEASE" ) {
   
    for(EERecHitCollection::const_iterator ite= eeHandle->begin(); ite != eeHandle->end(); ++ite ) { 
      
      EEDetId tmp_id(ite->id());
      float rh_calib = regionalCalibration_->getCalibMap()->coeff(tmp_id); 

      // calibrated energy
      float tmp_ene = ite->energy();
      float calib_ene = tmp_ene * rh_calib;

      // new rechit
      EcalRecHit aHit( *ite );
      aHit.setEnergy(calib_ene);
      EENewRecHits->push_back( aHit );
    }
  }

  // put the collection of reconstructed hits with previous round calibration in the event      
  LogInfo("CalibRechitsForZProducer") << "total # EB calibrated rechits: " << EBNewRecHits->size();
  LogInfo("CalibRechitsForZProducer") << "total # EE calibrated rechits: " << EENewRecHits->size();
  LogInfo("CalibRechitsForZProducer") << "total # EB rechits before calibration: " << ebHandle->size();
  LogInfo("CalibRechitsForZProducer") << "total # EE rechits before calibration: " << eeHandle->size();

  iEvent.put( EBNewRecHits, EBNewRecHitCollection_ );
  iEvent.put( EENewRecHits, EENewRecHitCollection_ );
}

//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE( CalibRechitsForZProducer );
