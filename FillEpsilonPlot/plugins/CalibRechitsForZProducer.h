#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibCode/CalibTools/interface/EcalRegionalCalibration.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

using namespace reco;

class CalibRechitsForZProducer : public edm::EDProducer {
   public:
      explicit CalibRechitsForZProducer(const edm::ParameterSet&);
      ~CalibRechitsForZProducer();
      virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);

   private:
      edm::Handle< EBRecHitCollection > ebHandle;
      edm::Handle< EBRecHitCollection > eeHandle;

      edm::EDGetTokenT<EBRecHitCollection> EBRecHitCollectionToken_;
      edm::EDGetTokenT<EERecHitCollection> EERecHitCollectionToken_;

      std::string EBNewRecHitCollection_;
      std::string EENewRecHitCollection_;

      EcalRegionalCalibrationBase *regionalCalibration_;
      EcalRegionalCalibration<EcalCalibType::Xtal> xtalCalib;

      int currentIteration_;
      std::string calibMapPath_; 
      bool isCRAB_;
      std::string Barrel_orEndcap_;
};
