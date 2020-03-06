import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.inputFiles = 'root://maite.iihe.ac.be//store/user/tomc/heavyNeutrinoMiniAOD/Moriond17_aug2018_miniAODv3/displaced/HeavyNeutrino_lljj_M-10_V-0.00107238052948_e_Dirac_massiveAndCKM_LO/heavyNeutrino_10.root'
options.parseArguments()

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag (process.GlobalTag, 'auto:run2_mc')

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
)

## Output file
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.OUT = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('test.root'),
    outputCommands = cms.untracked.vstring(['keep *'])
)

from PhysicsTools.PatAlgos.tools.jetTools import *

updateJetCollection(
    process,
    jetSource=cms.InputTag("slimmedJets"),
    jetCorrections=("AK4PFchs", cms.vstring(["L1FastJet", "L2Relative", "L3Absolute"]), "None"),
    btagInfos = ['pfImpactParameterTagInfos', 'pfInclusiveSecondaryVertexFinderTagInfos', 'pfDeepCSVTagInfos']
)

process.updatedPatJets.addBTagInfo = cms.bool(True)
process.updatedPatJets.addDiscriminators = cms.bool(True)
process.updatedPatJets.addJetCorrFactors = cms.bool(True)
process.updatedPatJets.addTagInfos = cms.bool(True)


process.globalFeatures = cms.EDProducer("JetFeatureXTagInfoProducer")
process.cpfFeatures = cms.EDProducer("ChargedCandidateXTagInfoProducer")
process.npfFeatures = cms.EDProducer("NeutralCandidateXTagInfoProducer")
process.svFeatures = cms.EDProducer("SecondaryVertexXTagInfoProducer")



process.task = cms.Task(
        process.pfImpactParameterTagInfos,
        process.pfInclusiveSecondaryVertexFinderTagInfos,
        process.pfDeepCSVTagInfos,
        process.patJetCorrFactors,
        process.globalFeatures,
        process.cpfFeatures,
        process.npfFeatures,
        process.svFeatures,
        process.updatedPatJets
)
process.p = cms.Path(process.task)

process.endpath= cms.EndPath(process.OUT)
