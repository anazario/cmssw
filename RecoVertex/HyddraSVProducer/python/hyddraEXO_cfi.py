import FWCore.ParameterSet.Config as cms

hyddraSVsEXOProducer = cms.EDProducer("HyddraSVsEXOProducer",
    muons        = cms.InputTag("slimmedMuons"),
    electrons    = cms.InputTag("slimmedElectrons"),
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    leptonic = cms.PSet(
        seedCosThetaCut = cms.double(-1.0),
        maxNormChi2     = cms.double(5.0),
    ),
)
