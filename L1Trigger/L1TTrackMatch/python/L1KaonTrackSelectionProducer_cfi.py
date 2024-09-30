import FWCore.ParameterSet.Config as cms

L1KaonTrackSelectionProducer = cms.EDProducer('L1KaonTrackSelectionProducer',
  #l1TracksInputTag = cms.InputTag("L1TrackNullSelectionProducer","Level1TTTracksNullSelected"),
  l1TracksInputTag = cms.InputTag("l1tTrackSelectionProducer","Level1TTTracksSelected"),                                            
  outputCollectionName = cms.string("Level1TTKaonTracksSelected"),
  processSimulatedTracks = cms.bool(True), # return selected tracks after cutting on the floating point values
  processEmulatedTracks = cms.bool(True), # return selected tracks after cutting on the bitwise emulated values
                                              debug = cms.int32(2) # Verbosity levels: 0, 1, 2, 3, 4
)

L1KaonTrackSelectionProducerExtended = L1KaonTrackSelectionProducer.clone(
  l1TracksInputTag = cms.InputTag("L1TrackSelectionProducerExtended", "Level1TTTracksExtendedSelected"),
  #l1TracksInputTag = cms.InputTag("L1TrackNullSelectionProducerExtended", "Level1TTTracksExtendedNullSelected"),
  outputCollectionName = "Level1TTKaonTracksExtendedSelected",
)

