#ifndef L1Trigger_L1TTrackMatch_L1TrackWordUnpacker_HH
#define L1Trigger_L1TTrackMatch_L1TrackWordUnpacker_HH

#include <vector>
#include <string>

#include <ap_fixed.h>
#include <ap_int.h>

#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

namespace l1trackunpacker {
  using L1TTTrackType = TTTrack<Ref_Phase2TrackerDigi_>;
  
  // For precision studies
  const unsigned int PT_INTPART_BITS{9}; 
  const unsigned int ETA_INTPART_BITS{3};
  const unsigned int kExtraGlobalPhiBit{4};

  using pt_intern     = ap_ufixed<TTTrack_TrackWord::TrackBitWidths::kRinvSize - 1, PT_INTPART_BITS, AP_TRN, AP_SAT>;
  using glbeta_intern = ap_fixed<TTTrack_TrackWord::TrackBitWidths::kTanlSize, ETA_INTPART_BITS, AP_TRN, AP_SAT>;
  using glbphi_intern = ap_int<TTTrack_TrackWord::TrackBitWidths::kPhiSize + kExtraGlobalPhiBit>;
  using z0_intern     = ap_int<TTTrack_TrackWord::TrackBitWidths::kZ0Size>;  // 40cm / 0.1
  using d0_intern     = ap_uint<TTTrack_TrackWord::TrackBitWidths::kD0Size>;

  const unsigned int DoubleToBit(double value, unsigned int maxBits, double step);
  const double BitToDouble(unsigned int bits, unsigned int maxBits, double step);
  double FloatPtFromBits(const L1TTTrackType& track);
  double FloatEtaFromBits(const L1TTTrackType& track);
  double FloatPhiFromBits(const L1TTTrackType& track);
  double FloatZ0FromBits(const L1TTTrackType& track);
}  // namespace l1trackunpacker

#endif
