import FWCore.ParameterSet.Config as cms

hyddraLeptonTracks = cms.EDProducer("HyddraLeptonTrackProducer",
    leptonType = cms.string("muon"),
    src        = cms.InputTag("slimmedMuons"),
)

hyddraSVsEXOProducer = cms.EDProducer("HyddraSVsEXOProducer",
    tracks       = cms.InputTag("hyddraLeptonTracks"),
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    leptonic = cms.PSet(
        seedCosThetaCut     = cms.double(-1.0),
        applySeedChi2Cut    = cms.bool(False),
        maxNormChi2         = cms.double(5.0),
        applyDcaCut         = cms.bool(False),
        maxDca              = cms.double(15.0),
        useSmoothing        = cms.bool(True),
        useMuonSystemBounds = cms.bool(True),
    ),
)
