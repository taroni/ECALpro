import FWCore.ParameterSet.Config as cms

ecalRecalRecHit = cms.EDProducer('CalibRechitsForZProducer',

    # parameters relevant for the calibration                                      
    CurrentIteration = cms.untracked.int32(0),                                      
    calibMapPath = cms.untracked.string('dummyMap.root'),
    isCRAB = cms.untracked.bool(False),
    Barrel_orEndcap = cms.untracked.string('Barrel_or_Endcap'),

    # input rechit collections
    EBRecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEB"),
    EERecHitCollectionTag = cms.untracked.InputTag("ecalRecHit","EcalRecHitsEE"),
 
    # output rechit collections
    EBNewRecHitCollection = cms.string('EcalRecalRecHitsEB'),
    EENewRecHitCollection = cms.string('EcalRecalRecHitsEE')
)
