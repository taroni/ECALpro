#ifndef CALIBELECTRONFORZ_H
#define CALIBELECTRONFORZ_H

#include <TROOT.h>
#include <TLorentzVector.h>

#include <vector>
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"


namespace calib
{
  class CalibElectronForZ {
    
  public:
    
    CalibElectronForZ();

  CalibElectronForZ(const reco::GsfElectron* ele,const EcalRecHitCollection* theHits, const EcalRecHitCollection* theEEHits) : 
    theElectron_(ele),
      theHits_(theHits), 
      theEEHits_(theEEHits) 
      {
      };

  CalibElectronForZ(const reco::GsfElectron* ele,const reco::SuperCluster* psc,const EcalRecHitCollection* theHits, const EcalRecHitCollection* theEEHits) : 
    theElectron_(ele),
      theParentSC_(psc),
      theHits_(theHits), 
      theEEHits_(theEEHits) 
      {
      };
    
    ~CalibElectronForZ() {};
    
    std::vector< std::pair<int,float> > getCalibModulesWeights(TString calibtype);
    const reco::GsfElectron* getRecoElectron() { return theElectron_; }
    const reco::SuperCluster* getParentSuperCluster() { return theParentSC_; }
    const EcalRecHitCollection* getRecHits() { return theHits_; }
    const EcalRecHitCollection* getEERecHits() { return theEEHits_; }
    
  private:
    
    const reco::GsfElectron* theElectron_;
    const reco::SuperCluster* theParentSC_;
    
    const EcalRecHitCollection* theHits_;
    const EcalRecHitCollection* theEEHits_;

  };
}
#endif

