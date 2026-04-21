import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('maxNormChi2',
                 5.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Max vertex chi2/ndof (default: 5.0)")
options.setDefault('maxEvents', -1)
options.setDefault('outputFile', 'hyddraEXO_output.root')
options.parseArguments()

process = cms.Process("HYDDRAEXO")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles or ['file:input.root']))

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

process.load("RecoVertex.HyddraSVProducer.hyddraEXO_cfi")
process.hyddraSVsEXOProducer.leptonic.maxNormChi2 = cms.double(options.maxNormChi2)

process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=cms.untracked.vstring(
        'drop *',
        'keep recoVertexs_hyddraSVsEXOProducer_seedVertices_*',
        'keep recoVertexs_hyddraSVsEXOProducer_inclusiveVertices_*',
        'keep recoVertexs_hyddraSVsEXOProducer_isolatedVertices_*',
        'keep ints_hyddraSVsEXOProducer_disambiguationFlags_*',
        'keep ints_hyddraSVsEXOProducer_seedIsolationFlags_*',
        'keep ints_hyddraSVsEXOProducer_isolationFlags_*',
        'keep recoTracks_hyddraSVsEXOProducer_leptonTracks_*',
    ),
)

process.p = cms.Path(process.hyddraSVsEXOProducer)
process.ep = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.p, process.ep)
