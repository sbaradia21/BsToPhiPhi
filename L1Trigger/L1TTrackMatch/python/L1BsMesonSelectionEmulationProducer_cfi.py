import FWCore.ParameterSet.Config as cms

L1BsMesonSelectionEmulationProducer = cms.EDProducer('L1BsMesonSelectionEmulationProducer',
  l1PhiMesonWordInputTag = cms.InputTag("L1PhiMesonSelectionEmulationProducer","Level1TTPhiMesonSelectedEmulation"),
  posTrackInputTag = cms.InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedEmulationPositivecharge"),
  negTrackInputTag = cms.InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedEmulationNegativecharge"),
  GenParticleInputTag = cms.InputTag("genParticles",""),
  outputCollectionName = cms.string("Level1TTBsMesonSelectedEmulation"),
  cutSet = cms.PSet(
           dzmax = cms.double(0.5),
          dxymax = cms.double(999.0),
           dRmin = cms.double(0.2),
           dRmax = cms.double(1.0),
    trkPairdRmax = cms.double(0.2),
     phiPairMmin = cms.double(5.0),
     phiPairMmax = cms.double(5.8),
         bsPtMin = cms.double(13.0)
  ),
  isSignal = cms.bool(True),
  debug = cms.int32(0)
)
