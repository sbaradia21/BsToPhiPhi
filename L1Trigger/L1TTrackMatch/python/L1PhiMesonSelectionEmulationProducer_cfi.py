import FWCore.ParameterSet.Config as cms
#from L1Trigger.VertexFinder.VertexProducer_cff import VertexProducer

L1PhiMesonSelectionEmulationProducer = cms.EDProducer('L1PhiMesonSelectionEmulationProducer',
  l1PosKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationPositivecharge"),
  l1NegKaonTracksInputTag = cms.InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationNegativecharge"),
  #L1VertexInputTag = cms.InputTag("VertexProducer", VertexProducer.l1VertexCollectionName.value()),                   
  outputCollectionName = cms.string("Level1TTPhiMesonSelectedEmulation"),
  cutSet = cms.PSet(
        dxymax = cms.double(999.0),
         dzmax = cms.double(0.5),
         dRmax = cms.double(0.5),
    tkpairMmin = cms.double(1.0),
    tkpairMmax = cms.double(1.03)
  ),
  debug = cms.int32(0)
)
