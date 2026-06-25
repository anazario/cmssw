import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('leptonType',
                 'muon',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Lepton type: 'muon' or 'electron'")
options.register('applySeedChi2Cut',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Apply maxNormChi2 to seed vertices")
options.register('maxNormChi2',
                 5.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Max vertex chi2/ndof (default: 5.0)")
options.register('applyDcaCut',
                 False,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Reject seeds with successful DCA calculation above maxDca")
options.register('maxDca',
                 15.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Maximum two-track DCA in cm when applyDcaCut is true")
options.register('useSmoothing',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Use Kalman vertex smoothing/refitted tracks")
options.register('useMuonSystemBounds',
                 True,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Extend Kalman vertex fitter bounds to the muon system")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraEXO_output.root')
options.parseArguments()

process = cms.Process("HYDDRAEXO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles or ['file:input.root']))

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

process.load("RecoVertex.HyddraSVProducer.hyddraEXO_cfi")

_src = 'slimmedMuons' if options.leptonType == 'muon' else 'slimmedElectrons'
process.hyddraLeptonTracks.leptonType = cms.string(options.leptonType)
process.hyddraLeptonTracks.src        = cms.InputTag(_src)
process.hyddraSVsEXOProducer.leptonic.applySeedChi2Cut = cms.bool(options.applySeedChi2Cut)
process.hyddraSVsEXOProducer.leptonic.maxNormChi2 = cms.double(options.maxNormChi2)
process.hyddraSVsEXOProducer.leptonic.applyDcaCut = cms.bool(options.applyDcaCut)
process.hyddraSVsEXOProducer.leptonic.maxDca = cms.double(options.maxDca)
process.hyddraSVsEXOProducer.leptonic.useSmoothing = cms.bool(options.useSmoothing)
process.hyddraSVsEXOProducer.leptonic.useMuonSystemBounds = cms.bool(options.useMuonSystemBounds)

process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep recoTracks_hyddraLeptonTracks__*',
        'keep recoVertexs_hyddraSVsEXOProducer_seedVertices_*',
        'keep recoVertexs_hyddraSVsEXOProducer_inclusiveVertices_*',
        'keep recoVertexs_hyddraSVsEXOProducer_isolatedVertices_*',
        'keep ints_hyddraSVsEXOProducer_disambiguationFlags_*',
        'keep ints_hyddraSVsEXOProducer_seedIsolationFlags_*',
        'keep ints_hyddraSVsEXOProducer_isolationFlags_*',
    ),
)

process.p = cms.Path(process.hyddraLeptonTracks + process.hyddraSVsEXOProducer)
process.ep = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p, process.ep)
