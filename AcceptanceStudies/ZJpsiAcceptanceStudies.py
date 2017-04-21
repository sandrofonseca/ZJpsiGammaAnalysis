import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/PhaseIFall16MiniAOD/ZToJPsiGamma-TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/FlatPU28to62HcalNZSRAW_PhaseIFall16_90X_upgrade2017_realistic_v6_C1-v2/80000/4E05E30E-B814-E711-B02F-0025905A607E.root'
    )
)

process.TFileService = cms.Service ('TFileService',
    fileName = cms.string ('acceptance.root')
)



process.acceptance = cms.EDAnalyzer('AcceptanceStudies',
    verbose = cms.bool(True),
    configName = cms.string("ZJpsiGamma 2017 MC sample"),	
    minMuPt = cms.double(0.0),# in GeV
    maxMuEta = cms.double(2.5),
    minMuonLeadPt = cms.double(5.0),# in GeV
    minMuonTrailPt = cms.double(2.0), # in GeV
    GammaMinPtCut = cms.double(12.0),# in GeV
    DeltaRLeadMuPhotonSel = cms.double(1.0),# deltaR>DeltaRLeadMuPhotonSel
    DeltaRTrailPhotonSel  = cms.double(1.0),# deltaR>DeltaRTrailPhotonSel 
    minJPsiMass = cms.double(2.95),# in GeV
    maxJPsiMass = cms.double(3.25)# in GeV

)


process.p = cms.Path(process.acceptance)
