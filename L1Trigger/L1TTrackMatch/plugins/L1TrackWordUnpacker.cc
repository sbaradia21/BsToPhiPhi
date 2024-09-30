#include <iostream>
#include <cmath>
#include <string>

#include "L1TrackWordUnpacker.h"

namespace l1trackunpacker {
  const unsigned int DoubleToBit(double value, unsigned int maxBits, double step) {
    unsigned int digitized_value = std::floor(std::abs(value) / step);
    unsigned int digitized_maximum = (1 << (maxBits - 1)) - 1;  // The remove 1 bit from nBits to account for the sign
    if (digitized_value > digitized_maximum)
      digitized_value = digitized_maximum;
    if (value < 0)
      digitized_value = (1 << maxBits) - digitized_value;  // two's complement encoding
    return digitized_value;
  }
  const double BitToDouble(unsigned int bits, unsigned int maxBits, double step) {
    int isign = 1;
    unsigned int digitized_maximum = (1 << maxBits) - 1;
    if (bits & (1 << (maxBits - 1))) {  // check the sign
      isign = -1;
      bits = (1 << (maxBits + 1)) - bits;
    }
    return (double(bits & digitized_maximum) + 0.5) * step * isign;
  }
  double FloatPtFromBits(const L1TTTrackType& track) {
    ap_uint<TTTrack_TrackWord::TrackBitWidths::kRinvSize - 1> ptBits = track.getRinvWord();
    pt_intern digipt;
    digipt.V = ptBits.range();
    return digipt.to_double();
  }
  double FloatEtaFromBits(const L1TTTrackType& track) {
    TTTrack_TrackWord::tanl_t etaBits = track.getTanlWord();
    glbeta_intern digieta;
    digieta.V = etaBits.range();
    return digieta.to_double();
  }
  double FloatPhiFromBits(const L1TTTrackType& track) {
    int sector = track.phiSector();
    double sector_phi_value = 0;
    if (sector < 5) {
      sector_phi_value = 2.0 * M_PI * sector / 9.0;
    } else {
      sector_phi_value = (-1.0 * M_PI + M_PI / 9.0 + (sector - 5) * 2.0 * M_PI / 9.0);
    }
    glbphi_intern trkphiSector = DoubleToBit(sector_phi_value, TTTrack_TrackWord::TrackBitWidths::kPhiSize + kExtraGlobalPhiBit, TTTrack_TrackWord::stepPhi0);
    glbphi_intern local_phiBits = 0;
    local_phiBits.V = track.getPhiWord();

    glbphi_intern local_phi =
      DoubleToBit(BitToDouble(local_phiBits, TTTrack_TrackWord::TrackBitWidths::kPhiSize, TTTrack_TrackWord::stepPhi0),
		  TTTrack_TrackWord::TrackBitWidths::kPhiSize + l1trackunpacker::kExtraGlobalPhiBit,
		  TTTrack_TrackWord::stepPhi0);
    glbphi_intern digiphi = local_phi + trkphiSector;
    return BitToDouble(digiphi, TTTrack_TrackWord::TrackBitWidths::kPhiSize + l1trackunpacker::kExtraGlobalPhiBit, TTTrack_TrackWord::stepPhi0);
  }
  double FloatZ0FromBits(const L1TTTrackType& track) {
    z0_intern trkZ = track.getZ0Word();
    return BitToDouble(trkZ, TTTrack_TrackWord::TrackBitWidths::kZ0Size, TTTrack_TrackWord::stepZ0);
  }
}
