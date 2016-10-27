#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibCode/CalibTools/interface/EcalRegionalCalibration.h"
#include "CalibCode/CalibTools/interface/CalibElectronForZ.h"
#include "CalibCode/FillEpsilonPlot/interface/JSON.h"
#include "Calibration/EcalCalibAlgos/interface/ZeeKinematicTools.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;
using namespace edm;
using namespace reco;

enum calibGranularity{ xtal, tt, etaring };

class FillEpsilonPlotForZ : public edm::EDAnalyzer {
   public:
      explicit FillEpsilonPlotForZ(const edm::ParameterSet&);
      ~FillEpsilonPlotForZ();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      
      static inline float invMassCalc(float Energy1, float Eta1, float Phi1, float Energy2, float Eta2, float Phi2) {
	return (sqrt(2 * Energy1 * Energy2 * (1 - cosTheta12(Eta1, Phi1, Eta2, Phi2))));
      }

      static inline float cosTheta12(float Eta1, float Phi1, float Eta2, float Phi2) {
	return ((cos(Phi1 - Phi2) + sinh(Eta1) * sinh(Eta2)) / (cosh(Eta1) * cosh(Eta2)));
      }

      static const int N_XTAL_TOTAL  = 75848;
      static const int N_XTAL_BARREL = 61200;
      static const int N_XTAL_ENDCAP = 14648;

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ---------- user defined ------------------------
      // void getWeight(float recomass, calib::CalibElectronForZ* ele, float evweight, float mcWeight);
      void getWeight(float recomass, calib::CalibElectronForZ* ele, float evweight, float mcWeight, int subDetId, const EcalRecHitCollection*);
      void bookHistograms();
      TH1F** initializeEpsilonHistograms(const char *name, const char *title, int size, int isweight );  
      void deleteEpsilonPlot(TH1F **h, int size); 
      void writeEpsilonPlot(TH1F **h, int size);        

      // ----------member data ---------------------------
      bool requireOppositeCharge_;
      double minInvMassCut_, maxInvMassCut_;
      int electronSelection_;
      double elePtCut_, eleEtaCut_;
      double maxDReleSc_;
      bool useMassInsteadOfEpsilon_;
      std::string massMethod_;

      bool isMC_;
      bool isCRAB_;
      std::string Barrel_orEndcap_;
      std::string outputDir_;
      std::string outfilename_;
      std::string calibMapPath_; 
      int currentIteration_;
      double SystOrNot_;
      std::string jsonFile_; 

      // int channels_;
      int nRegionsEB_, nRegionsEE_;

      calibGranularity calibTypeNumber_;
      EcalRegionalCalibration<EcalCalibType::Xtal> xtalCalib;
      EcalRegionalCalibrationBase *regionalCalibration_;

      // ---------- collections ---------------------------
      edm::EDGetTokenT<EBRecHitCollection> EBRecHitCollectionToken_;
      edm::EDGetTokenT<EERecHitCollection> EERecHitCollectionToken_;
      edm::EDGetTokenT<SuperClusterCollection> EBSuperClusterCollectionToken_;
      edm::EDGetTokenT<SuperClusterCollection> EESuperClusterCollectionToken_;
      edm::EDGetTokenT<GsfElectronCollection> ElectronCollectionToken_;
      std::string mcProducer_;

      edm::Handle< EBRecHitCollection > hits;
      edm::Handle< EERecHitCollection > ehits;
      edm::Handle< SuperClusterCollection > ebScCollection;
      edm::Handle< SuperClusterCollection > eeScCollection;
      edm::Handle< GsfElectronCollection > electronCollection;


      // ---------- outputs ---------------------------
      TFile *outfile_;

      TH1F **weightedRescaleFactorEB;
      TH1F **unweightedRescaleFactorEB;
      TH1F **weightEB;
      TH1F **weightedRescaleFactorEE;
      TH1F **unweightedRescaleFactorEE;
      TH1F **weightEE;
      
      TH2F *h2_xtalRecalibCoeffEB;
      TH2F *h2_xtalRecalibCoeffEEp;
      TH2F *h2_xtalRecalibCoeffEEm;
      TH1F *h1_residualMiscalibEB;
      TH1F *h1_residualMiscalibEEp;
      TH1F *h1_residualMiscalibEEm;
      TH1F *EventFlow;
      TH1F *allEpsilon_EE; 
      TH1F *allEpsilon_EEnw; 
      TH1F *allEpsilon_EB;
      TH1F *allEpsilon_EBnw;
      TH2F *entries_EEp;
      TH2F *entries_EEm;
      TH2F *entries_EB;
      TH2F *Occupancy_EEp;
      TH2F *Occupancy_EEm;
      TH2F *Occupancy_EB;
      TH2F *zMassVsIetaEB;
      TH2F *zMassVsETEB;

      TTree *myTree;
      int   isEBEB;
      int   isEEEE;
      int   isHR9HR9;
      float mass4tree;


      float thisEventW;

      // Json file
      JSON* myjson;
};
