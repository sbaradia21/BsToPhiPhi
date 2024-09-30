// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1PhiMesonSelectionEmulationProducer
//
/**\class L1PhiMesonSelectionEmulationProducer L1PhiMesonSelectionEmulationProducer.cc L1Trigger/L1TTrackMatch/plugins/L1PhiMesonSelectionEmulationProducer.cc

 Description: Build Phi meson candidates from positively and negatively charged selected L1Tracks (kaons)

 Implementation:
     Inputs:
         std::vector<TTTrack> - Each floating point TTTrack inside this collection inherits from
                                a bit-accurate TTTrack_TrackWord, used for emulation purposes.
     Outputs:
         l1t::TkLightMesonWordCollection - A collection of phi meson candidates
*/
// ----------------------------------------------------------------------------                                                    
// Original Authors: Alexx Perloff, Pritam Palit                                                                                   
//         Created:  Thu, 16 Dec 2021 19:02:50 GMT                                                                                 
//          Authors: Sweta Baradia, Suchandra Dutta, Subir Sarkar (April 2024)                                                     
//----------------------------------------------------------------------------- 

// system include files
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <TMath.h>
#include <cmath>
#include <bitset>
#include <TSystem.h>
#include <array>
// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>
#include <stdio.h>
#include <cassert>
#include <cstdlib>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"
#include "DataFormats/L1Trigger/interface/TkLightMesonWord.h"
#include "L1TrackWordUnpacker.h"


using namespace std;
using namespace edm;
using namespace l1t;

using L1TTTrackType            = TTTrack<Ref_Phase2TrackerDigi_>;                                 
using TTTrackCollection        = std::vector<L1TTTrackType>;                                    
using TTTrackRef               = edm::Ref<TTTrackCollection>;                                   
using TTTrackRefCollection     = edm::RefVector<TTTrackCollection>;                             
using TTTrackCollectionHandle  = edm::Handle<TTTrackRefCollection>;
using TTTrackRefCollectionUPtr = std::unique_ptr<TTTrackRefCollection>; 

// Relevant constants for the converted track word    
enum TrackBitWidths {
  kPtSize = TTTrack_TrackWord::TrackBitWidths::kRinvSize -1,  // Width of pt
  kPtMagSize = 9,//9                                              // Width of pt magnitude (unsigned)
  kEtaSize = TTTrack_TrackWord::TrackBitWidths::kTanlSize,     // Width of eta //by pp
  kEtaMagSize = 3,                                             // Width of eta magnitude (signed) //3 by PP 
};


namespace l1tphimesonemu {
  const unsigned int kInternalPhiWidth{8};
  const unsigned int kGlobalPhiExtra{4};
  static constexpr double minPhi0{-0.7853981696};
  using global_phi_t = ap_uint<TTTrack_TrackWord::TrackBitWidths::kPhiSize>;
  //typedef ap_uint<15> global_phislice_t;
  const unsigned int kGlobalPhiBins = 1 << kInternalPhiWidth;
  //  const unsigned int kGlobalPhiTotalBins = 1 << TTTrack_TrackWord::TrackBitWidths::kPhiSize;
  //const unsigned int kGlobalPhiTotalBins = 1 << 12;
  const double kStepPhi = (2 * -minPhi0) / kGlobalPhiBins;
  const unsigned int kNSector{9};
  const unsigned int kNQuadrants{4};

  double undigitizeSignedValue(unsigned int twosValue, unsigned int nBits) { //Used when local to global phi converstion done using PP's way
    // Check that none of the bits above the nBits-1 bit, in a range of [0, nBits-1], are set.
    // This makes sure that it isn't possible for the value represented by `twosValue` to be
    // any bigger than ((1 << nBits) - 1).
    assert((twosValue >> nBits) == 0);
    
    // Convert from twos compliment to C++ signed integer (normal digitized value)
    int digitizedValue = twosValue;
    if (twosValue & (1 << (nBits - 1))) {  // check if the twosValue is negative
      digitizedValue -= (1 << nBits);
    }
    // Convert to floating point value
    return (double(digitizedValue) + 0.5);
  }

  void printTrackInfo(const L1TTTrackType& track, double phiEmu) {
    cout << "("
	 << track.momentum().perp() << ", " 
	 << track.momentum().eta() << ", " 
	 << track.momentum().phi() << ", "
	 << track.getStubRefs().size() << ", " 
	 << track.stubPtConsistency() << ", " 
	 << track.chi2ZRed() << ", "
	 << track.chi2XYRed() << ", " 
	 << track.z0() 
	 << ")" << endl;
   
    
    //Pt
    ap_uint<TrackBitWidths::kPtSize> ptEmulationBits = track.getRinvWord();
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_TRN, AP_SAT> ptEmulation;
    //ap_ufixed<10,7, AP_TRN, AP_SAT> ptEmulation;
    //    ptEmulation.V = ptEmulationBits.range(11,2);
    ptEmulation.V = ptEmulationBits.range();
    // Eta
    TTTrack_TrackWord::tanl_t etaEmulationBits = track.getTanlWord();
    ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize, AP_TRN, AP_SAT> etaEmulation;
    etaEmulation.V = etaEmulationBits.range();

    cout << "(" 
	<< ptEmulation.to_double() << ", " 
	<< etaEmulation.to_double() << ", " 
	<< track.getPhi() << ", "
        << track.getNStubs() << ", " 
	<< track.getBendChi2() << ", " 
	<< track.getChi2RZ() << ", " 
	<< track.getChi2RPhi() << ", " 
	<< track.getZ0() 
	<< ")" << endl;
  }
}


//Class declaration begins...
class L1PhiMesonSelectionEmulationProducer : public edm::global::EDProducer<> {
public:
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  explicit L1PhiMesonSelectionEmulationProducer(const edm::ParameterSet&);
  ~L1PhiMesonSelectionEmulationProducer() override;
  void beginJob() override;
  void endJob();

  static constexpr double KaonMass = 0.493677 ;
  static constexpr double ETAPHI_LSB = M_PI / (1 << 12);
  static constexpr double Z0_LSB = 0.05;
  void trackPropertiesBitVsFloat(const TTTrackCollectionHandle posTrackHandle, const TTTrackCollectionHandle negTrackHandle) const;

  //PP Local to Global Phi conversion ; not using now
  l1tphimesonemu::global_phi_t localToGlobalPhi(const TTTrack_TrackWord::phi_t& local_phi,
						const l1tphimesonemu::global_phi_t& sector_shift) const;
  
  std::vector<l1tphimesonemu::global_phi_t> generatePhiSliceLUT (unsigned int N) {  //Used when local to global phi converstion done using PP's way
    //    std::cout<<"Doing for N : "<<N<<std::endl;
    //std::cout<<"-------------------------"<<std::endl;
    float sliceCentre = 0.0;
    std::vector<l1tphimesonemu::global_phi_t> phiLUT;
    for (unsigned int q = 0; q <= N; q++) {
      phiLUT.push_back((l1tphimesonemu::global_phi_t)(sliceCentre / l1tphimesonemu::kStepPhi));
      //std::cout<<"q="<<q<<" sliceCentre is :"<<sliceCentre<<" and LUT val is : "<<phiLUT[q]<<" an kstepPhi was : "<<l1tphimesonemu::kStepPhi<<std::endl;
      sliceCentre += 2 * M_PI / N;
    }
    return phiLUT;
  }  
private:
  // ----------constants, enums and typedefs ---------
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
 
  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors
  // ----------member data ---------------------------
  const edm::EDGetTokenT<TTTrackRefCollection> posTrackToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> negTrackToken_;
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  const double dxymax_, dzmax_, dRmax_, tkpairMmin_, tkpairMmax_;

  int debug_;
  std::vector<l1tphimesonemu::global_phi_t> phiQuadrants; //Used when local to global phi converstion done using PP's way
  std::vector<l1tphimesonemu::global_phi_t> phiShift; //Used when local to global phi converstion done using PP's way

  mutable std::array<unsigned int, 8> evCounters_;
  mutable std::array<unsigned int, 5> evTrkCounters_;
  TH1D *evtCutFlow, *evtCutFlowTrk;
  TH1D *nTrkPos, *nTrkNeg, *nPhi;
  TH1D *dPtTrkPairH_, *dEtaTrkPairH_, *dPhiTrkPairH_; //Sim - Emu for each track
  TH1D *dz0TrkPairH_, *dRTrkPairH_, *dRTrkPair1H_, *dRTrkPair2H_,*dz0TrkPairSimH_, *dRTrkPairSimH_,*diff_dRTrkPair, *diff_dz0TrkPair; //Pos - Neg for Emu and Sim
  TH1D *poskptSimH_, *posketaSimH_, *poskphiSimH_, *negkptSimH_, *negketaSimH_, *negkphiSimH_, *poskptH_, *posketaH_, *poskphiH_, *negkptH_, *negketaH_, *negkphiH_; //quantities in emu and sim
  TH2D *poskpt2DH_, *negkpt2DH_, *diff_poskpt2D, *diff_negkpt2D;
  TH1D *diff_poskpt, *diff_posketa, *diff_poskphi, *diff_negkpt, *diff_negketa, *diff_negkphi;// *diff_poskSinHeta, *diff_negkSinHeta; //Sim -Emu for pos and neg tracks
  TH1D *diff_poskpx, *diff_poskpy, *diff_poskpz, *diff_poskE, *diff_negkpx, *diff_negkpy, *diff_negkpz, *diff_negkE; //Sim -Emu for pos and neg tracks
  TH1D *diff_kpx, *diff_kpy, *diff_kpz, *diff_kE; //Sim -Emu for pos and neg tracks
  TH1D *phimassH_,*phimassSimH_,*phimass1H_,*phimass1SimH_,*phimass2H_,*phimass2SimH_,*phimass3H_,*phimass3SimH_,*phimass4H_,*phimass4SimH_; // For Sim and Emu at dfrnt levels
  TH1D *phimass1_diff;
  
};

//
// constructors and destructor
//
L1PhiMesonSelectionEmulationProducer::L1PhiMesonSelectionEmulationProducer(const edm::ParameterSet& iConfig)
  : posTrackToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1PosKaonTracksInputTag"))),
    negTrackToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1NegKaonTracksInputTag"))),
    outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
    cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
    dxymax_(cutSet_.getParameter<double>("dxymax")),
    dzmax_(cutSet_.getParameter<double>("dzmax")),
    dRmax_(cutSet_.getParameter<double>("dRmax")),
    tkpairMmin_(cutSet_.getParameter<double>("tkpairMmin")),
    tkpairMmax_(cutSet_.getParameter<double>("tkpairMmax")),
    debug_(iConfig.getParameter<int>("debug")) {
  produces<l1t::TkLightMesonWordCollection>(outputCollectionName_);
    cout << "--> L1PhiMesonSelectionEmulationProducer:"
       << " dxymax: " << dxymax_
       << " dzmax: " << dzmax_
       << " dRmax: " << dRmax_
       << " tkpairMmin: " << tkpairMmin_
       << " tkpairMmax: " << tkpairMmax_
       << endl;
  phiQuadrants = L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(l1tphimesonemu::kNQuadrants); //Used when local to global phi converstion done using PP's way 
  phiShift = L1PhiMesonSelectionEmulationProducer::generatePhiSliceLUT(l1tphimesonemu::kNSector); //Used when local to global phi converstion done using PP's way 
}

L1PhiMesonSelectionEmulationProducer::~L1PhiMesonSelectionEmulationProducer() {}
void L1PhiMesonSelectionEmulationProducer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return;
  evtCutFlow   = fs->make<TH1D>("evtCutFlow", "",8,-0.5,7.5);
  for (size_t i = 0; i < evCounters_.size(); ++i) evCounters_[i] = 0;
  evtCutFlowTrk   = fs->make<TH1D>("evtCutFlowTrk", "",5,-0.5,4.5);
  for (size_t i = 0; i < evTrkCounters_.size(); ++i) evTrkCounters_[i] = 0;
  nTrkPos    = fs->make<TH1D>("nTrkPos", "Number of positive Tracks", 150, -0.5, 149.5);
  nTrkNeg    = fs->make<TH1D>("nTrkNeg", "Number of negative Tracks", 150, -0.5, 149.5);
  nPhi      = fs->make<TH1D>("nPhi", "Number of Phi meson candidates", 20, -0.5, 19.5);

  dPtTrkPairH_  = fs->make<TH1D>("dPtTrkPair", "dPt between Track pairs", 200, -0.1, 0.1);
  dEtaTrkPairH_ = fs->make<TH1D>("dEtaTrkPair", "dEta between Track pairs", 200, -0.05, 0.05);
  dPhiTrkPairH_ = fs->make<TH1D>("dPhiTrkPair", "dPhi between Track pairs", 200, -0.001, 0.001);
  
  dz0TrkPairH_  = fs->make<TH1D>("dz0TrkPair", "dz0 between Track pairs", 200, 0, 20);
  dRTrkPairH_   = fs->make<TH1D>("dRTrkPair", "dR between Track pairs", 200, 0.0, 10.0);
  dRTrkPair1H_   = fs->make<TH1D>("dRTrkPair1", "dR between Track pairs alternate formula", 200, 0.0, 10.0);
  dRTrkPair2H_  = fs->make<TH1D>("dRTrkPair2", "dR between Track pairs after Phi selection", 200, 0.0, 0.5);
  dz0TrkPairSimH_  = fs->make<TH1D>("dz0TrkPairSim", "dz0 between Track pairs Sim", 200, 0, 20);
  dRTrkPairSimH_   = fs->make<TH1D>("dRTrkPairSim", "dR between Track pairs Sim", 200, 0.0, 10.0);
  diff_dRTrkPair   = fs->make<TH1D>("diff_dRTrkPair", "dR between Track pairs Sim - Emu", 200, 0.0, 10.0);
  diff_dz0TrkPair  = fs->make<TH1D>("diff_dz0TrkPair", "dz0 between Track pairs Sim - Emu", 200, 0.0, 10.0);

  poskptSimH_ = fs->make<TH1D>("poskpt Sim", "pt",  320, 0, 10);
  posketaSimH_ = fs->make<TH1D>("posketa Sim", "eta",  100, -3, 3);
  poskphiSimH_ = fs->make<TH1D>("poskphi Sim", "phi",  100, -4, 4);
  negkptSimH_ = fs->make<TH1D>("negkpt Sim", "pt", 320, 0, 10);
  negketaSimH_ = fs->make<TH1D>("negketa Sim", "eta", 100, -3, 3);
  negkphiSimH_ = fs->make<TH1D>("negkphi Sim", "phi",  100, -4, 4);

  poskptH_ = fs->make<TH1D>("poskpt", "pt",  320, 0, 10);
  posketaH_= fs->make<TH1D>("posketa", "eta",  100, -3, 3);
  poskphiH_ = fs->make<TH1D>("poskphi", "phi",  100, -4, 4);
  negkptH_ = fs->make<TH1D>("negkpt", "pt", 320, 0, 10);
  negketaH_ = fs->make<TH1D>("negketa ", "eta", 100, -3, 3);
  negkphiH_ = fs->make<TH1D>("negkphi ", "phi",  100, -4, 4);

  poskpt2DH_ = fs->make<TH2D>("poskpt2D", "pt sim vs emu",  2000, 0, 100,2000,0,100);
  negkpt2DH_ = fs->make<TH2D>("negkpt2D", "pt sim vs emu",  2000, 0, 100,2000,0,100);
  
  //Sim - Emu
  diff_poskpt = fs->make<TH1D>("diff_poskpt", "Difference of pt", 200, -0.5, 0.5);
  diff_posketa = fs->make<TH1D>("diff_posketa", "Difference of eta", 200, -0.1, 0.1);
  //diff_poskSinHeta = fs->make<TH1D>("diff_poskSinHeta", "Difference of eta", 200, -0.1, 0.1);
  diff_poskphi = fs->make<TH1D>("diff_poskphi", "Difference of phi", 200, -0.05, 0.05);
  diff_negkpt = fs->make<TH1D>("diff_negkpt", "Difference of pt", 200, -0.5, 0.5);
  diff_negketa = fs->make<TH1D>("diff_negketa", "Difference of eta", 200, -0.1, 0.1);
  //diff_negkSinHeta = fs->make<TH1D>("diff_negkSinHeta", "Difference of eta", 200, -0.1, 0.1);
  diff_negkphi = fs->make<TH1D>("diff_negkphi", "Difference of phi", 200, -0.05, 0.05);

  diff_poskpt2D = fs->make<TH2D>("diff_poskpt2D", "Difference of pt vs sim pt", 200, -0.5, 0.5,2000,0,100);
  diff_negkpt2D = fs->make<TH2D>("diff_negkpt2D", "Difference of pt vs sim pt", 200, -0.5, 0.5,2000,0,100);

  diff_poskpx = fs->make<TH1D>("diff_poskpx", "Difference of px", 100, -2,2);
  diff_poskpy = fs->make<TH1D>("diff_poskpy", "Difference of py", 100, -2,2);
  diff_poskpz = fs->make<TH1D>("diff_poskpz", "Difference of pz", 100, -2,2);
  diff_poskE = fs->make<TH1D>("diff_poskE", "Difference of E", 100, -2,2);
  diff_negkpx = fs->make<TH1D>("diff_negkpx", "Difference of px", 100, -2,2);
  diff_negkpy = fs->make<TH1D>("diff_negkpy", "Difference of py", 100, -2,2);
  diff_negkpz = fs->make<TH1D>("diff_negkpz", "Difference of pz", 100, -2,2);
  diff_negkE = fs->make<TH1D>("diff_negkE", "Difference of E", 100, -2,2);
  diff_kpx = fs->make<TH1D>("diff_kpx", "Difference of px", 100, -2,2);
  diff_kpy = fs->make<TH1D>("diff_kpy", "Difference of py", 100, -2,2);
  diff_kpz = fs->make<TH1D>("diff_kpz", "Difference of pz", 100, -2,2);
  diff_kE = fs->make<TH1D>("diff_kE", "Difference of E", 100, -2,2);

  phimassH_ = fs->make<TH1D>("phimass","Final from Phiword", 1000, 0.98, 1.1);

  phimass1H_ = fs->make<TH1D>("phimass1","Before any cut", 1000, 0.98, 1.1);
  phimass1SimH_ = fs->make<TH1D>("phimass1Sim","Before any cut", 1000, 0.98, 1.1);
  phimass1_diff = fs->make<TH1D>("phimass1_diff","Before any cut", 1000, -1, 1);

  phimass2H_ = fs->make<TH1D>("phimass2","After dz and dR cut", 1000, 0.98, 1.1);
  phimass2SimH_ = fs->make<TH1D>("phimass2Sim","After dz and dR cut", 1000, 0.98, 1.1);

  phimass3H_ = fs->make<TH1D>("phimass3","After all cuts", 1000, 0.98, 1.1);
  phimass3SimH_ = fs->make<TH1D>("phimass3Sim","After all cuts", 1000, 0.98, 1.1);

  phimass4H_ = fs->make<TH1D>("phimass4","After all cuts bin",200, 0.98, 1.1);
  phimass4SimH_ = fs->make<TH1D>("phimass4Sim","After all cuts bin",200, 0.98, 1.1);
}

l1tphimesonemu::global_phi_t L1PhiMesonSelectionEmulationProducer::localToGlobalPhi(const TTTrack_TrackWord::phi_t& local_phi, const l1tphimesonemu::global_phi_t& sector_shift) const { //Used when local to global phi converstion done using PP's way 
  int PhiMin = 0;
  int PhiMax = phiQuadrants.back();
  
  int phiMultiplier = TTTrack_TrackWord::TrackBitWidths::kPhiSize - l1tphimesonemu::kInternalPhiWidth;
  
  int tempPhi = floor(l1tphimesonemu::undigitizeSignedValue(local_phi, TTTrack_TrackWord::TrackBitWidths::kPhiSize) / pow(2, phiMultiplier)) + sector_shift;
  //int tempPhi = floor(l1tphimesonemu::undigitizeSignedValue(local_phi, TTTrack_TrackWord::TrackBitWidths::kPhiSize)) + sector_shift;    
  //  int tempPhi = floor(l1tphimesonemu::undigitizeSignedValue(local_phi.to_uint(), TTTrack_TrackWord::TrackBitWidths::kPhiSize)) + sector_shift;  
  
  if (tempPhi < PhiMin) {
    tempPhi += PhiMax;
  } 
  else if (tempPhi > PhiMax) {
    tempPhi -= PhiMax;
  }                                                    
  return l1tphimesonemu::global_phi_t(tempPhi);
}

// ------------ method called to produce the data  ------------                                                                                                    
void L1PhiMesonSelectionEmulationProducer::trackPropertiesBitVsFloat(const TTTrackCollectionHandle posTrackHandle,
                                                                     const TTTrackCollectionHandle negTrackHandle) const {
  unsigned int nTracks = posTrackHandle->size() + negTrackHandle->size();
  for (size_t i = 0; i < nTracks; ++i) {
    const auto& trackRef = (i < posTrackHandle->size()) ? posTrackHandle->at(i) : negTrackHandle->at(i - posTrackHandle->size());
    const auto& track = *trackRef;
    dPtTrkPairH_->Fill(track.momentum().perp() - l1trackunpacker::FloatPtFromBits(track));
    dEtaTrkPairH_->Fill(track.momentum().eta() - l1trackunpacker::FloatEtaFromBits(track));
    dPhiTrkPairH_->Fill(track.momentum().phi() - l1trackunpacker::FloatPhiFromBits(track));
  }
}
                                                                                                                                                                     
void L1PhiMesonSelectionEmulationProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  auto L1PhiMesonEmulationOutput = std::make_unique<l1t::TkLightMesonWordCollection>();

  TTTrackCollectionHandle posTrackHandle;
  iEvent.getByToken(posTrackToken_, posTrackHandle);

  TTTrackCollectionHandle negTrackHandle;
  iEvent.getByToken(negTrackToken_, negTrackHandle);

  auto nPosKaonOutputApproximate = posTrackHandle->size();
  auto nNegKaonOutputApproximate = negTrackHandle->size();
  auto nPhiMesonOutputApproximate = nPosKaonOutputApproximate + nNegKaonOutputApproximate;

  L1PhiMesonEmulationOutput->reserve(nPhiMesonOutputApproximate);

  
  ap_ufixed<64, 32> etaphi_conv = 1.0 / ETAPHI_LSB;
  ap_ufixed<64, 32> z0_conv = 1.0 * Z0_LSB;

#if 0 
  cout << "==> L1PhiMesonSelectionEmulationProducer::produce" << endl;
#endif

  ++evCounters_[0];
  evtCutFlow->Fill(0);
  
  if(nPosKaonOutputApproximate > 1 && nNegKaonOutputApproximate > 1) {
    ++evCounters_[1];
    evtCutFlow->Fill(1);
  }
  nTrkPos->Fill(nPosKaonOutputApproximate);
  nTrkNeg->Fill(nNegKaonOutputApproximate);
  
  trackPropertiesBitVsFloat(posTrackHandle, negTrackHandle);
  std::array<unsigned int, 4> iCounters = {{0, 0, 0, 0}};
  
  for (size_t i = 0; i < nPosKaonOutputApproximate; i++) {
    const auto& trackPosKaonRef = posTrackHandle->at(i);
    const auto& trackPosKaon = *trackPosKaonRef;
    
    //Using Unpacker
    //double trkptPos  = l1trackunpacker::FloatPtFromBits(trackPosKaon);
    //double trketaPos = l1trackunpacker::FloatEtaFromBits(trackPosKaon);
    //double trkphiPos = l1trackunpacker::FloatPhiFromBits(trackPosKaon);
    //double trkz0Pos  = l1trackunpacker::FloatZ0FromBits(trackPosKaon);

    ap_uint<TrackBitWidths::kPtSize> ptEmulationBitsPos = trackPosKaon.getRinvWord();
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_TRN, AP_SAT> ptEmulationPos;
    //ap_ufixed<10,7,AP_TRN, AP_SAT> ptEmulationPos;
    //    ptEmulationPos.V = ptEmulationBitsPos.range(11,2);
    ptEmulationPos.V = ptEmulationBitsPos.range();
    double trkptPos = ptEmulationPos.to_double();//+0.015625;//+(2^-3)/2
    //double trkptPos = trackPosKaon.momentum().perp(); //Taking pt from simulation

   TTTrack_TrackWord::tanl_t etaEmulationBitsPos = trackPosKaon.getTanlWord();
   ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize, AP_TRN, AP_SAT> etaEmulationPos;
   etaEmulationPos.V = etaEmulationBitsPos.range();
   double trketaPos = etaEmulationPos.to_double();

   l1tphimesonemu::global_phi_t trkphiEmuPos = L1PhiMesonSelectionEmulationProducer::localToGlobalPhi(trackPosKaon.getPhiWord(), phiShift[trackPosKaon.phiSector()]);
   //   double trkphiPos = trkphiEmuPos * l1tphimesonemu::kStepPhi;
   //   if (trkphiPos > M_PI) trkphiPos -= 2 * M_PI;
   double trkphiPos = l1trackunpacker::FloatPhiFromBits(trackPosKaon); //Taking Phi from Triplet Unpacker

#if 0    
    cout << "1. phi: " << trackPosKaon.momentum().phi() << ", phiemu: " << trkphiPos << endl; 
#endif

    double trkz0Pos = trackPosKaon.getZ0();
    double trkpxPos = trkptPos * cos(trkphiPos);
    double trkpyPos = trkptPos * sin(trkphiPos);
    double trkpzPos = trkptPos * sinh(trketaPos);
    double trkePos  = sqrt(pow(trkptPos,2) + pow(trkpzPos,2) + pow(KaonMass,2));

#if 0
    cout << "+ve track" << endl;
    l1tphimesonemu::printTrackInfo(trackPosKaon);
#endif
    
    // Global  variables from track objects for positive tracks
    float l1postkpt = trackPosKaon.momentum().perp(); 
    float l1postketa = trackPosKaon.momentum().eta();
    float l1postkphi = trackPosKaon.momentum().phi();
    float l1postkpx = l1postkpt*cos(l1postkphi);
    float l1postkpy = l1postkpt*sin(l1postkphi);
    float l1postkpz = l1postkpt*sinh(l1postketa);
    float l1postke  = sqrt(pow(l1postkpt,2) + pow(l1postkpz,2) + pow(KaonMass,2));
    //SB GLobal variables from track word for positive tracks
    float l1emupostkpt =  trkptPos;
    float l1emupostketa = trketaPos;
    float l1emupostkphi = (trkphiPos > M_PI) ? trkphiPos - 2 * M_PI : trkphiPos;
    float diff_poskptG = l1postkpt - l1emupostkpt;
    float diff_posketaG = l1postketa - l1emupostketa;
    float diff_poskphiG = l1postkphi - l1emupostkphi;
    //float diff_poskSinHetaG = sinh(l1postketa) - sinh(l1emupostketa);
        
    poskptSimH_->Fill(l1postkpt);
    posketaSimH_->Fill(l1postketa);
    poskphiSimH_->Fill(l1postkphi);
    poskptH_->Fill(l1emupostkpt);
    poskpt2DH_->Fill(l1postkpt,l1emupostkpt);
    posketaH_->Fill(l1emupostketa);
    poskphiH_->Fill(l1emupostkphi);
    diff_poskpt->Fill(diff_poskptG);
    diff_poskpt2D->Fill(diff_poskptG,l1postkpt);
    diff_posketa->Fill(diff_posketaG);
    //    diff_poskSinHeta->Fill(diff_poskSinHetaG);
    diff_poskphi->Fill(diff_poskphiG);
    diff_poskpx->Fill(trkpxPos-l1postkpx);
    diff_poskpy->Fill(trkpyPos-l1postkpy);
    diff_poskpz->Fill(trkpzPos-l1postkpz);
    diff_poskE->Fill(trkePos-l1postke);
    
    for (size_t j = 0; j < nNegKaonOutputApproximate; j++) {
      ++evTrkCounters_[0];
      evtCutFlowTrk->Fill(0);
      const auto& trackNegKaonRef = negTrackHandle->at(j);
      const auto& trackNegKaon = *trackNegKaonRef;
      
      //Using Unpacker
      //double trkptNeg  = l1trackunpacker::FloatPtFromBits(trackNegKaon);
      //double trketaNeg = l1trackunpacker::FloatEtaFromBits(trackNegKaon);
      //double trkphiNeg = l1trackunpacker::FloatPhiFromBits(trackNegKaon);
      //double trkz0Neg = l1trackunpacker::FloatZ0FromBits(trackNegKaon);

      ap_uint<TrackBitWidths::kPtSize> ptEmulationBitsNeg = trackNegKaon.getRinvWord();
      ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize, AP_TRN, AP_SAT> ptEmulationNeg;
      //ap_ufixed<10,7, AP_TRN, AP_SAT> ptEmulationNeg;
      //ptEmulationNeg.V = ptEmulationBitsNeg.range(11,2);
      ptEmulationNeg.V = ptEmulationBitsNeg.range();
      double trkptNeg = ptEmulationNeg.to_double();//+0.015625;
      ///double trkptNeg = trackNegKaon.momentum().perp(); //Taking pt from simulation 
      
      TTTrack_TrackWord::tanl_t etaEmulationBitsNeg = trackNegKaon.getTanlWord();
      ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize, AP_TRN, AP_SAT> etaEmulationNeg;
      etaEmulationNeg.V = etaEmulationBitsNeg.range();
      double trketaNeg = etaEmulationNeg.to_double();

      l1tphimesonemu::global_phi_t trkphiEmuNeg = L1PhiMesonSelectionEmulationProducer::localToGlobalPhi(trackNegKaon.getPhiWord(), phiShift[trackNegKaon.phiSector()]);
      //double trkphiNeg = trkphiEmuNeg * l1tphimesonemu::kStepPhi;
      //if (trkphiNeg > M_PI) trkphiNeg -= 2 * M_PI;
      double trkphiNeg =  l1trackunpacker::FloatPhiFromBits(trackNegKaon);
      
     #if 0
      cout << "2. phi: " << trackNegKaon.momentum().phi() << ", phiemu: " << trkphiNeg << endl; 
#endif

      double trkz0Neg = trackNegKaon.getZ0();      
      double trkpxNeg = trkptNeg * cos(trkphiNeg);
      double trkpyNeg = trkptNeg * sin(trkphiNeg);
      double trkpzNeg = trkptNeg * sinh(trketaNeg);
      double trkeNeg  = sqrt(pow(trkptNeg,2) + pow(trkpzNeg,2) + pow(KaonMass,2));

#if 0
      cout << "-ve track" << endl;
      l1tphimesonemu::printTrackInfo(trackNegKaon);
#endif

      
      //SB Global  variables of track objects for negetive tracks                                                                                
      float l1negtkpt = trackNegKaon.momentum().perp(); 
      float l1negtketa = trackNegKaon.momentum().eta();
      float l1negtkphi = trackNegKaon.momentum().phi();
      float l1negtkpx = l1negtkpt*cos(l1negtkphi);
      float l1negtkpy = l1negtkpt*sin(l1negtkphi);
      float l1negtkpz = l1negtkpt*sinh(l1negtketa);
      float l1negtke  = sqrt(pow(l1negtkpt,2) + pow(l1negtkpz,2) + pow(KaonMass,2));
      
      //SB GLobal phi from track word  for negetive tracks                                                                                            
      float l1emunegtkpt = trkptNeg;
      float l1emunegtketa = trketaNeg;
      float l1emunegtkphi = (trkphiNeg > M_PI) ? trkphiNeg - 2 * M_PI : trkphiNeg;
      float diff_negkptG = l1negtkpt - l1emunegtkpt;
      float diff_negketaG = l1negtketa - l1emunegtketa;
      float diff_negkphiG = l1negtkphi - l1emunegtkphi;
      //float diff_negkSinHetaG = sinh(l1negtketa) - sinh(l1emunegtketa);
      
      if(i == 0) {
	negkptSimH_->Fill(l1negtkpt);
	negketaSimH_->Fill(l1negtketa);
	negkphiSimH_->Fill(l1negtkphi);
	negkptH_->Fill(l1emunegtkpt);
	negkpt2DH_->Fill(l1negtkpt,l1emunegtkpt);
	negketaH_->Fill(l1emunegtketa);
	negkphiH_->Fill(l1emunegtkphi);
	diff_negkpt->Fill(diff_negkptG);
	diff_negkpt2D->Fill(diff_negkptG,l1negtkpt);
	diff_negketa->Fill(diff_negketaG);
	//diff_negkSinHeta->Fill(diff_negkSinHetaG);
	diff_negkphi->Fill(diff_negkphiG);
	diff_negkpx->Fill(trkpxNeg-l1negtkpx);
	diff_negkpy->Fill(trkpyNeg-l1negtkpy);
	diff_negkpz->Fill(trkpzNeg-l1negtkpz);
	diff_negkE->Fill(trkeNeg-l1negtke);
      }

      //Getting sum of px, py, pz and e for both tracks using track object variables 
      float lpxPhi = l1negtkpx + l1postkpx;
      float lpyPhi = l1negtkpy + l1postkpy;
      float lpzPhi = l1negtkpz + l1postkpz;
     float lePhi = l1negtke + l1postke;
      //Getting sum of px, py, pz and e for both tracks using trackword info
      double pxPhi = trkpxNeg + trkpxPos;
      double pyPhi = trkpyNeg + trkpyPos;
      double pzPhi = trkpzNeg + trkpzPos;
      double ePhi  = trkeNeg  + trkePos;
      
      diff_kpx->Fill(pxPhi-lpxPhi);
      diff_kpy->Fill(pyPhi-lpyPhi);
      diff_kpz->Fill(pzPhi-lpzPhi);
      diff_kE->Fill(ePhi-lePhi);

      double trkmasspairPhi = sqrt(pow(ePhi,2) - pow(pxPhi, 2) - pow(pyPhi, 2) - pow(pzPhi, 2)); //Emulation                                                    
      double phimassPair = sqrt(pow(lePhi,2) - pow(lpxPhi, 2) - pow(lpyPhi, 2) - pow(lpzPhi, 2)); //Simulation 
      
      //math::PtEtaPhiMLorentzVectorD kaon1_vec(trkptPos, trketaPos, trkphiPos, KaonMass);
      //math::PtEtaPhiMLorentzVectorD kaon2_vec(trkptNeg, trketaNeg, trkphiNeg, KaonMass);
      //double trkmasspairPhiVec  = (kaon1_vec+kaon2_vec).M();
      
      phimass1SimH_->Fill(phimassPair);  //Before any cut
      phimass1H_->Fill(trkmasspairPhi);
      phimass1_diff->Fill(trkmasspairPhi-phimassPair);
      
      //dzo
      dz0TrkPairSimH_->Fill(std::fabs(trackPosKaon.POCA().z() - trackNegKaon.POCA().z()));
      dz0TrkPairH_->Fill(std::fabs(trkz0Pos - trkz0Neg));
      if (std::fabs(trkz0Pos - trkz0Neg) > dzmax_) continue;
      ++iCounters[0];
      ++evTrkCounters_[1];
      evtCutFlowTrk->Fill(1);

      //DR cut
      double convdPhi = std::abs(trkphiPos - trkphiNeg);
      if(convdPhi > M_PI) convdPhi-= 2 * M_PI;
      double trkdrpairPhi = sqrt(pow(convdPhi,2) + pow(trketaPos - trketaNeg, 2));
      double trkdrpairPhi1 = reco::deltaR(trketaPos, trkphiPos, trketaNeg, trkphiNeg);
      dRTrkPairH_->Fill(trkdrpairPhi); //DR_emu
      dRTrkPair1H_->Fill(trkdrpairPhi1);

      math::PtEtaPhiMLorentzVector itrkposP4(trackPosKaon.momentum().perp(), trackPosKaon.momentum().eta(), trackPosKaon.momentum().phi(), KaonMass);
      math::PtEtaPhiMLorentzVector itrknrgP4(trackNegKaon.momentum().perp(), trackNegKaon.momentum().eta(), trackNegKaon.momentum().phi(), KaonMass);
      dRTrkPairSimH_->Fill(reco::deltaR(itrkposP4,itrknrgP4)); //DR_sim

      if (trkdrpairPhi > dRmax_) continue;
      ++iCounters[1];
      ++evTrkCounters_[2];
      evtCutFlowTrk->Fill(2);

      diff_dz0TrkPair->Fill(std::fabs(std::fabs(trkz0Pos - trkz0Neg) - std::fabs(trackPosKaon.POCA().z() - trackNegKaon.POCA().z())));
      diff_dRTrkPair->Fill(std::fabs(trkdrpairPhi-reco::deltaR(itrkposP4,itrknrgP4)));

      phimass2SimH_->Fill(phimassPair); //After dz0 and DR cut
      phimass2H_->Fill(trkmasspairPhi);

      //mass cut
      if (trkmasspairPhi < tkpairMmin_ || trkmasspairPhi > tkpairMmax_) continue;
      ++iCounters[2];
      ++evTrkCounters_[3];
      evtCutFlowTrk->Fill(3);
      
      l1t::TkLightMesonWord::valid_t validPhi           = trackPosKaon.getValid() && trackNegKaon.getValid();
      l1t::TkLightMesonWord::pt_t ptPhi                 = sqrt(pow(pxPhi, 2) + pow(pyPhi, 2)); 
      l1t::TkLightMesonWord::glbphi_t phiPhi            = atan2(pyPhi, pxPhi) / ETAPHI_LSB;
      l1t::TkLightMesonWord::glbeta_t etaPhi            = asinh(pzPhi/sqrt(pow(pxPhi, 2) + pow(pyPhi, 2))) / ETAPHI_LSB;
      l1t::TkLightMesonWord::z0_t z0Phi                 = ((trkz0Pos + trkz0Neg) / Z0_LSB)* 0.5;
      l1t::TkLightMesonWord::mass_t mPhi                = sqrt(pow(ePhi, 2) - pow(pxPhi, 2) - pow(pyPhi, 2) - pow(pzPhi, 2));
      l1t::TkLightMesonWord::type_t typePhi             = l1t::TkLightMesonWord::TkLightMesonTypes::kPhiType;
      l1t::TkLightMesonWord::ntracks_t ntracksPhi       = 2;
      l1t::TkLightMesonWord::index_t firstTrkIndex      = i;
      l1t::TkLightMesonWord::index_t secondTrkIndex     = j;
      l1t::TkLightMesonWord::unassigned_t unassignedPhi = 0;

#if 0
      cout << "1. phi: " << trackPosKaon.momentum().phi() << ", phiemu: " << ((trkphiPos > M_PI) ? trkphiPos - 2 * M_PI : trkphiPos) << endl; 
      cout << "2. phi: " << trackNegKaon.momentum().phi() << ", phiemu: " << ((trkphiNeg > M_PI) ? trkphiNeg - 2 * M_PI : trkphiNeg) << endl; 
      cout << "Phi Candidate write: ("
	   << ptPhi << ", " 
	   << etaPhi << ", " 
	   << phiPhi << ", "
	   << mPhi
	   << ")" << endl;
#endif      
      l1t::TkLightMesonWord PhiWord(validPhi, ptPhi, phiPhi, etaPhi, z0Phi, mPhi, typePhi, ntracksPhi, firstTrkIndex, secondTrkIndex, unassignedPhi);
      bool dupl = false;
      for (const auto& el: *L1PhiMesonEmulationOutput) {
        double ptDiff  = el.pt()  - PhiWord.pt();
	double etaDiff = el.glbeta() - PhiWord.glbeta();
        double phiDiff = el.glbphi() - PhiWord.glbphi();
        if ( fabs(etaDiff) < 1.0e-03 &&
             fabs(phiDiff) < 1.0e-03 &&
             fabs(ptDiff)  < 1.0e-02 )
          {
            dupl = true;
#if 0
            cout << "--> Duplicate Phi found!" << endl;
#endif
            break;
          }
      }
      if(!dupl) {
	++iCounters[3];
	++evTrkCounters_[4];
	evtCutFlowTrk->Fill(4);
	phimass3SimH_->Fill(phimassPair); //After all cuts
	phimass3H_->Fill(trkmasspairPhi);
	phimassH_->Fill(mPhi.to_double());
	phimass4SimH_->Fill(phimassPair); //After all cuts
        phimass4H_->Fill(trkmasspairPhi);

	dRTrkPair2H_->Fill(trkdrpairPhi);
	L1PhiMesonEmulationOutput->push_back(PhiWord);
      }
    }
  }

  for (size_t i = 0; i < iCounters.size(); ++i) {
    if (iCounters[i] > 0) {
      ++evCounters_[2+i];
      evtCutFlow->Fill(2+i);
    }
  }


  if(L1PhiMesonEmulationOutput->size() > 0) {
    ++evCounters_[6];
    evtCutFlow->Fill(6);
  }

  if(L1PhiMesonEmulationOutput->size() > 1) {
    ++evCounters_[7];
    evtCutFlow->Fill(7);
  }

  nPhi->Fill(L1PhiMesonEmulationOutput->size());
  iEvent.put(std::move(L1PhiMesonEmulationOutput), outputCollectionName_);

}

void L1PhiMesonSelectionEmulationProducer::endJob() {
  ostringstream dzTrkCutTag;
  dzTrkCutTag << "dzTrkPair < " << fixed << setprecision(1) << dzmax_ << " cm";

  ostringstream drTrkCutTag;
  drTrkCutTag << "drTrkPair < " << fixed << setprecision(1) << dRmax_ ;

  ostringstream dmTrkPairCutTag;
  dmTrkPairCutTag << fixed << setprecision(1) << tkpairMmin_ << " < Phi Mass < " << tkpairMmax_ << " GeV";

  std::vector<std::string> tags {
      "#events passing the filter",
      " #PosTrk > 1 & #NegTrk > 1",
      dzTrkCutTag.str(),
      drTrkCutTag.str(),
      dmTrkPairCutTag.str(),
      "Non Dupl Phi",
      "                  #phi > 0",
      "                  #phi > 1",
  };
  std::vector<std::string> tagsTrk {
      "#Trk Pairs",
      dzTrkCutTag.str(),
      drTrkCutTag.str(),
      dmTrkPairCutTag.str(),
      "Non Dupl Phi",
   };

  /*  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"----------Cutflow Phi-Evt Emulation----------"<<std::endl;
  std::cout<<"---------------------------------------------"<<std::endl;

  for (size_t i = 0; i < tags.size(); ++i)
    std::cout << setw(40) << tags[i] << setw(10) << evCounters_[i] << std::endl;
  /*  std::cout<<"----------Cutflow Phi-Trk Emulation----------"<<std::endl;
  for (size_t i = 0; i < tagsTrk.size(); ++i)
    std::cout << setw(40) << tagsTrk[i] << setw(10) << evTrkCounters_[i] << std::endl;
  */

  int nbins = evtCutFlow->GetNbinsX();
  for (int i = 1; i <=nbins; ++i) {
    evtCutFlow->GetXaxis()->SetBinLabel(i, tags[i-1].c_str());
  }

  int nbinsTrk = evtCutFlowTrk->GetNbinsX();
  for (int i = 1; i <=nbinsTrk; ++i) {
    evtCutFlowTrk->GetXaxis()->SetBinLabel(i, tagsTrk[i-1].c_str());
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1PhiMesonSelectionEmulationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PosKaonTracksInputTag",
                          edm::InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationPositivecharge"));
  desc.add<edm::InputTag>("l1NegKaonTracksInputTag",
                          edm::InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedEmulationNegativecharge"));
  desc.add<std::string>("outputCollectionName", "Level1TTPhiMesonSelectedEmulation");

  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dxymax", 1.0)->setComment("dxy must be less than this value, [cm]");
    descCutSet.add<double>("dzmax", 1.0)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("dRmax", 0.12)->setComment("dr must be less than this value, []");
    descCutSet.add<double>("tkpairMmin", 1.0)->setComment("tkpair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("tkpairMmax", 1.03)->setComment("tkpair mass must be less than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);
  }
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1PhiMesonSelectionEmulationProducer);
