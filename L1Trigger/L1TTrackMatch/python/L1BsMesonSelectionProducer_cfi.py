import FWCore.ParameterSet.Config as cms

L1BsMesonSelectionProducer = cms.EDProducer('L1BsMesonSelectionProducer',
  l1PhiCandInputTag = cms.InputTag("L1PhiMesonSelectionProducer", "Level1TTPhiMesonSelected"),
  outputCollectionName = cms.string("Level1TTBsMesonSelected"),
  cutSet = cms.PSet(
      dxymax = cms.double(999.0),
      dzmax = cms.double(0.5),
      dRmin = cms.double(0.2),
      dRmax = cms.double(1.0),
      trkpairdRmax = cms.double(0.2),
      phipairMmin = cms.double(5.0),
      phipairMmax = cms.double(5.8),
      ptmin = cms.double(13.0)
  ),
  isSignal = cms.bool(True),
  debug = cms.int32(0)
)

