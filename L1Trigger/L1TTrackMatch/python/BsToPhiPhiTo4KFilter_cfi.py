import FWCore.ParameterSet.Config as cms

BsToPhiPhiFilter = cms.EDFilter('BsToPhiPhiTo4KFilter',
  GenParticleInputTag = cms.InputTag("genParticles",""),
  filter = cms.bool(True)
)
