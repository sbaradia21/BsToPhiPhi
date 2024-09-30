import FWCore.ParameterSet.Config as cms

L1PhiMesonSelectionProducer = cms.EDProducer('L1PhiMesonSelectionProducer',
  l1PosKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedPositivecharge"),
  l1NegKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedNegativecharge"),
  outputCollectionName = cms.string("Level1TTPhiMesonSelected"),
  cutSet = cms.PSet(
      dxymax = cms.double(999.0),
      dzmax = cms.double(0.5),
      dRmax = cms.double(0.5),
      tkpairMmin = cms.double(1.0),
      tkpairMmax = cms.double(1.03)
  ),
  debug = cms.int32(0)
)


