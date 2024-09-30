// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1KaonTrackSelectionProducer
//
/**\class L1KaonTrackSelectionProducer L1KaonTrackSelectionProducer.cc L1Trigger/L1TTrackMatch/plugins/L1KaonTrackSelectionProducer.cc

 Description: Selects two set of positively and negatively charged L1Tracks corresponding to Kaons which already passed the criteria for for Light Meson track selection

 Implementation:
     Inputs:
         std::vector<TTTrack> - Each floating point TTTrack inside this collection inherits from
                                a bit-accurate TTTrack_TrackWord, used for emulation purposes.
     Outputs:
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack properties
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack_TrackWord properties
*/
//
// Original Author:  Alexx Perloff
//         Created:  Thu, 16 Dec 2021 19:02:50 GMT
//
//

// system include files
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "CommonTools/Utils/interface/AndSelector.h"
#include "CommonTools/Utils/interface/EtaRangeSelector.h"
#include "CommonTools/Utils/interface/MinSelector.h"
#include "CommonTools/Utils/interface/MinFunctionSelector.h"
#include "CommonTools/Utils/interface/MinNumberSelector.h"
#include "CommonTools/Utils/interface/PtMinSelector.h"
#include "CommonTools/Utils/interface/Selection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1D.h>
#include <TFile.h>
//
// class declaration
//
class L1KaonTrackSelectionProducer : public edm::global::EDProducer<> {
public:
  explicit L1KaonTrackSelectionProducer(const edm::ParameterSet&);
  ~L1KaonTrackSelectionProducer() override;
  void beginJob() override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word
  enum TrackBitWidths {
    kPtSize = TTTrack_TrackWord::TrackBitWidths::kRinvSize - 1,  // Width of pt
    kPtMagSize = 9,                                              // Width of pt magnitude (unsigned)
    kEtaSize = TTTrack_TrackWord::TrackBitWidths::kTanlSize,     // Width of eta
    kEtaMagSize = 3,                                             // Width of eta magnitude (signed)
  };

  typedef TTTrack<Ref_Phase2TrackerDigi_> L1Track;
  typedef std::vector<L1Track> TTTrackCollection;
  typedef edm::Ref<TTTrackCollection> TTTrackRef;
  typedef edm::RefVector<TTTrackCollection> TTTrackRefCollection;
  typedef edm::Handle<TTTrackRefCollection> TTTrackCollectionHandle;
  typedef std::unique_ptr<TTTrackRefCollection> TTTrackRefCollectionUPtr;

  // ----------member functions ----------------------
  void printDebugInfo(const TTTrackCollectionHandle& l1TracksHandle,
                      const TTTrackRefCollectionUPtr& vTTTrackOutput,
                      const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const;
  void printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation = false) const;
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors
  // Track charge selection 

  bool TTTrackChargeSelector(const L1Track& t) const { /*std::cout << "track r inv : " << t.rInv() << "    track charge : " << std::signbit(t.rInv()) << std::endl;*/ return std::signbit(t.rInv()) ;}; 

  bool TTTrackWordChargeSelector (const L1Track& t) const {
    //std::cout << "track wordd charge : " << std::signbit(t.getRinv()) << std::endl;
    //    return std::signbit(t.getRinv());

    ap_uint<1> chargeEmulationBits = t.getTrackWord()(
						      TTTrack_TrackWord::TrackBitLocations::kRinvMSB, TTTrack_TrackWord::TrackBitLocations::kRinvMSB);

    //edm::LogInfo("L1KaonTrackSelectionProducer") << "produce::Emulation track properties::ap_uint(bits) = " << ptEmulationBits.to_string(2)
    //                                         << " ap_ufixed(bits) = " << ptEmulation.to_string(2) << " ap_ufixed(float) = " << ptEmulation.to_double();
    //std::cout << "track word charge : " << chargeEmulationBits.to_uint() << std::endl;
    // You can also use directly std::signbit(getRinv) => getRinv gives you the signed real value of Rinv 
    return chargeEmulationBits.to_uint();

  };

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TTTrackRefCollection> l1TracksToken_;
  const std::string outputCollectionName_;
  bool processSimulatedTracks_, processEmulatedTracks_;
  int debug_;
  TH1D *z0TrkH_;
  TH1D *z0TrkPOCAH_;

};

//
// constructors and destructor
//
L1KaonTrackSelectionProducer::L1KaonTrackSelectionProducer(const edm::ParameterSet& iConfig)
  : l1TracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1TracksInputTag"))),
      outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
      processSimulatedTracks_(iConfig.getParameter<bool>("processSimulatedTracks")),
      processEmulatedTracks_(iConfig.getParameter<bool>("processEmulatedTracks")),
      debug_(iConfig.getParameter<int>("debug")) {
  // Confirm the the configuration makes sense
  if (!processSimulatedTracks_ && !processEmulatedTracks_) {
    throw cms::Exception("You must process at least one of the track collections (simulated or emulated).");
  }
 
  // Get additional input tags and define the EDM output based on the previous configuration parameters
  if (processSimulatedTracks_) {
    produces<TTTrackRefCollection>(outputCollectionName_ + "Positivecharge");
    produces<TTTrackRefCollection>(outputCollectionName_ + "Negativecharge");
  }
  if (processEmulatedTracks_) {
    produces<TTTrackRefCollection>(outputCollectionName_ + "EmulationPositivecharge");
    produces<TTTrackRefCollection>(outputCollectionName_ + "EmulationNegativecharge");
  }
}

L1KaonTrackSelectionProducer::~L1KaonTrackSelectionProducer() {}

//
// member functions
//
void L1KaonTrackSelectionProducer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return;
  z0TrkH_ = fs->make<TH1D>("z0Trk", "z0", 150, -15, 15);
  z0TrkPOCAH_ = fs->make<TH1D>("z0TrkPOCA", "POCA().z()", 150, -15, 15);
}
void L1KaonTrackSelectionProducer::printDebugInfo(const TTTrackCollectionHandle& l1TracksHandle,
						  const TTTrackRefCollectionUPtr& vTTTrackOutput,
						  const TTTrackRefCollectionUPtr& vTTTrackEmulationOutput) const {
  edm::LogInfo log("L1KaonTrackSelectionProducer");
  log << "The original track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are ... \n";
  for (const auto& track : *l1TracksHandle) {
    printTrackInfo(log, track, debug_ >= 4);
  }
  log << "\t---\n\tNumber of tracks in this selection = " << l1TracksHandle->size() << "\n\n";
  if (processSimulatedTracks_) {
    log << "The selected track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are ... \n";
    for (const auto& track : *vTTTrackOutput) {
      printTrackInfo(log, track, debug_ >= 4);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << vTTTrackOutput->size() << "\n\n";
  }
  if (processEmulatedTracks_) {
    log << "The emulation selected track collection (pt, eta, phi, nstub, bendchi2, chi2rz, chi2rphi, z0) values are "
           "... \n";
    for (const auto& track : *vTTTrackEmulationOutput) {
      printTrackInfo(log, track, debug_ >= 4);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << vTTTrackEmulationOutput->size() << "\n\n";
  }
  if (processSimulatedTracks_ && processEmulatedTracks_) {
    TTTrackRefCollection inSimButNotEmu;
    TTTrackRefCollection inEmuButNotSim;
    std::set_difference(vTTTrackOutput->begin(),
                        vTTTrackOutput->end(),
                        vTTTrackEmulationOutput->begin(),
                        vTTTrackEmulationOutput->end(),
                        std::back_inserter(inSimButNotEmu));
    std::set_difference(vTTTrackEmulationOutput->begin(),
                        vTTTrackEmulationOutput->end(),
                        vTTTrackOutput->begin(),
                        vTTTrackOutput->end(),
                        std::back_inserter(inEmuButNotSim));
    log << "The set of tracks selected via cuts on the simulated values which are not in the set of tracks selected "
           "by cutting on the emulated values ... \n";
    for (const auto& track : inSimButNotEmu) {
      printTrackInfo(log, track, debug_ >= 3);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << inSimButNotEmu.size() << "\n\n"
        << "The set of tracks selected via cuts on the emulated values which are not in the set of tracks selected "
           "by cutting on the simulated values ... \n";
    for (const auto& track : inEmuButNotSim) {
      printTrackInfo(log, track, debug_ >= 3);
    }
    log << "\t---\n\tNumber of tracks in this selection = " << inEmuButNotSim.size() << "\n\n";
  }
}
void L1KaonTrackSelectionProducer::printTrackInfo(edm::LogInfo& log, const TTTrackRef& track, bool printEmulation) const {
  log << "\t("
      << track->momentum().perp() << ", " 
      << track->momentum().eta() << ", " 
      << track->momentum().phi() << ", "
      << track->getStubRefs().size() << ", " 
      << track->stubPtConsistency() << ", " 
      << track->chi2ZRed() << ", "
      << track->chi2XYRed() << ", " 
      << track->z0() 
      << ")\n";

  if (printEmulation) {
    // pT
    ap_uint<TrackBitWidths::kPtSize> ptEmulationBits = track->getTrackWord()(
        TTTrack_TrackWord::TrackBitLocations::kRinvMSB - 1, TTTrack_TrackWord::TrackBitLocations::kRinvLSB);
    ap_ufixed<TrackBitWidths::kPtSize, TrackBitWidths::kPtMagSize> ptEmulation;
    ptEmulation.V = ptEmulationBits.range();
    // Eta
    TTTrack_TrackWord::tanl_t etaEmulationBits = track->getTanlWord();
    ap_fixed<TrackBitWidths::kEtaSize, TrackBitWidths::kEtaMagSize> etaEmulation;
    etaEmulation.V = etaEmulationBits.range();
    log << "\t\t(" 
	<< ptEmulation.to_double() << ", " 
	<< etaEmulation.to_double() << ", " 
	<< track->getPhi() << ", "
        << track->getNStubs() << ", " 
	<< track->getBendChi2() << ", " 
	<< track->getChi2RZ() << ", " 
	<< track->getChi2RPhi() << ", " 
	<< track->getZ0() 
	<< ")\n";
  }
}
// ------------ method called to produce the data  ------------
void L1KaonTrackSelectionProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  auto vTTPosTrackOutput = std::make_unique<TTTrackRefCollection>();
  auto vTTPosTrackEmulationOutput = std::make_unique<TTTrackRefCollection>();

  auto vTTNegTrackOutput = std::make_unique<TTTrackRefCollection>();
  auto vTTNegTrackEmulationOutput = std::make_unique<TTTrackRefCollection>();

  TTTrackCollectionHandle l1TracksHandle;

  iEvent.getByToken(l1TracksToken_, l1TracksHandle);
  size_t nOutputApproximate = l1TracksHandle->size();
  if (processSimulatedTracks_) {
    vTTPosTrackOutput->reserve(nOutputApproximate);
    vTTNegTrackOutput->reserve(nOutputApproximate);
  }
  if (processEmulatedTracks_) {
    vTTPosTrackEmulationOutput->reserve(nOutputApproximate);
    vTTNegTrackEmulationOutput->reserve(nOutputApproximate);
  }

  for (size_t i = 0; i < nOutputApproximate; i++) {
    const auto& trackRef = l1TracksHandle->at(i);
    const auto& track = *trackRef;
    z0TrkH_->Fill(trackRef->z0());
    z0TrkPOCAH_->Fill(trackRef->POCA().z());

    // Select tracks based on the floating point TTTrack
    if (processSimulatedTracks_) {
      (!TTTrackChargeSelector(track) ? vTTPosTrackOutput->push_back(trackRef)
             	                     : vTTNegTrackOutput->push_back(trackRef));
    }
    // Select tracks based on the bitwise accurate TTTrack_TrackWord
    if (processEmulatedTracks_) {
      (!TTTrackWordChargeSelector(track) ? vTTPosTrackEmulationOutput->push_back(trackRef)
                                         : vTTNegTrackEmulationOutput->push_back(trackRef));
    }
    /*  if (debug_ >= 2) {
	  printDebugInfo(l1TracksHandle,	vTTPosTrackOutput, TTPosTrackEmulationOutput);
	  printDebugInfo(l1TracksHandle,	vTTNegTrackOutput, vTTNegTrackEmulationOutput);
	}*/
  }  
  // Put the outputs into the event
  if (processSimulatedTracks_) {
    iEvent.put(std::move(vTTPosTrackOutput), outputCollectionName_ + "Positivecharge");
    iEvent.put(std::move(vTTNegTrackOutput), outputCollectionName_ + "Negativecharge");
  }
  if (processEmulatedTracks_) {
    iEvent.put(std::move(vTTPosTrackEmulationOutput), outputCollectionName_ + "EmulationPositivecharge");
    iEvent.put(std::move(vTTNegTrackEmulationOutput), outputCollectionName_ + "EmulationNegativecharge");
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1KaonTrackSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1KaonTrackSelectionProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1TracksInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<std::string>("outputCollectionName", "Level1TTKaonTracksSelected");
  desc.add<bool>("processSimulatedTracks", true)
      ->setComment("return selected tracks after cutting on the floating point values");
  desc.add<bool>("processEmulatedTracks", true)
      ->setComment("return selected tracks after cutting on the bitwise emulated values");
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1KaonTrackSelectionProducer);
