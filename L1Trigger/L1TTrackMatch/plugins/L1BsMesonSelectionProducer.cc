// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1BsMesonSelectionProducer
//
/**\class L1BsMesonSelectionProducer L1BsMesonSelectionProducer.cc L1Trigger/L1TTrackMatch/plugins/L1BsMesonSelectionProducer.cc

 Description: Selects two set of positively and negatively charged L1Tracks corresponding to Kaons which already passed the criteria for for Light Meson track selection

 Implementation:
     Inputs:
         std::vector<TTTrack> - Each floating point TTTrack inside this collection inherits from
                                a bit-accurate TTTrack_TrackWord, used for emulation purposes.
     Outputs:
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack properties
         std::vector<TTTrack> - A collection of TTTracks selected from cuts on the TTTrack_TrackWord properties
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
#include <array>

// user include files
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkBsCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidateFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkBsCandidateFwd.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Gen-level stuff
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include <TH1D.h>
#include <TFile.h>

//
// class declaration

//
using namespace std;
using namespace edm;
using namespace l1t;

class L1BsMesonSelectionProducer : public edm::one::EDProducer<edm::one::SharedResources> {
//class L1BsMesonSelectionProducer : public edm::global::EDProducer<> {
public:
  using L1TTTrackType            = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection        = std::vector<L1TTTrackType>;
  using TTTrackRefCollection     = edm::RefVector<TTTrackCollection>;
  using TTTrackCollectionHandle  = edm::Handle<TTTrackRefCollection>;
  using L1TTStubCollection =
    std::vector<edm::Ref<edmNew::DetSetVector<TTStub<Ref_Phase2TrackerDigi_>>, TTStub<Ref_Phase2TrackerDigi_>>>;

  explicit L1BsMesonSelectionProducer(const edm::ParameterSet&);
  ~L1BsMesonSelectionProducer() override;
  void beginJob() override;
  void endJob() override;
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static constexpr double KaonMass = 0.493677; // GeV

  static bool duplicateTrk(edm::Ptr<L1TTTrackType> trka, edm::Ptr<L1TTTrackType> trkb, const TrackerTopology& tTopo);
  static void printHeader(int ntrk);
  static void printTrk(edm::Ptr<L1TTTrackType> trkPtr, const TrackerTopology& tTopo, int itrk);
  static void stubInfo(edm::Ptr<L1TTTrackType> trk,
		       const TrackerTopology& tTopo,
		       int& nStub,
		       int& nStubSS,
		       int& nStubPS);
  static void fillTrackQuantities(edm::Ptr<L1TTTrackType> trk, 
				  const TrackerTopology& tTopo, 
				  TH1D* trkPtH, 
				  TH1D* trkEtaH, 
				  TH1D* trkPhiH, 
				  TH1D* trkZ0H, 
				  TH1D* trknStubH, 
				  TH1D* trknStubPSH,
				  TH1D* trkChi2BendH, 
				  TH1D* trkChi2ZH, 
				  TH1D* trkChi2XYH);
  bool isGenMatched(const reco::GenParticleCollection& genParticles,
                    const vector<TkPhiCandidate>& phiList,
                    bool verbose=false);

    
private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word

  using TkPhiCandidateCollectionHandle = edm::Handle<TkPhiCandidateCollection>;

  // ----------member functions ----------------------
  // void produce(edm::StreamID, edm::Event&, const edm::EventSetup&)  const override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TkPhiCandidateCollection> l1PhiCandToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> posTrackToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> negTrackToken_;

  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  const std::string outputCollectionName_;
  const bool isSignal_;
  int debug_;
  const edm::ParameterSet cutSet_;
  const double dxymax_, dzmax_, dRmin_, dRmax_, trkPairdRmax_, phiPairMmin_, phiPairMmax_, bsPtMin_;

  TH1D *ptAllPhiH_, *etaAllPhiH_, *phiAllPhiH_, *massAllPhiH_, *z0AllPhiH_;
  TH1D *nPhiH_, *nBsH_,*dzPhiPairH_, *drPhiPairH_, *drPhi1TrkPairH_, *drPhi2TrkPairH_, *bsMassH_, *bsMassFH_, *bsPtH_, *bsEtaH_, *bsPhiH_;
  TH1D *bsMassBinH_, *bsMassFBinH_;
  mutable std::array<unsigned int, 11> evCounters_;
  TH1D *evtCutFlow_;
  std::array<TH1D*, 4> trkPtH_, trkEtaH_, trkPhiH_, trkZ0H_, trknStubH_, trknStubPSH_, trkChi2BendH_, trkChi2ZH_, trkChi2XYH_;
  std::array<TH1D*, 2> ptPhiH_, etaPhiH_, phiPhiH_, massPhiH_, z0PhiH_;
  TH1D *drGenPhiH_, *dptGenPhiH_, *detaGenPhiH_;
};

//
// constructors and destructor
//
L1BsMesonSelectionProducer::L1BsMesonSelectionProducer(const edm::ParameterSet& iConfig)
    : l1PhiCandToken_(consumes<TkPhiCandidateCollection>(iConfig.getParameter<edm::InputTag>("l1PhiCandInputTag"))),
    tTopoToken_(esConsumes<TrackerTopology, TrackerTopologyRcd>(edm::ESInputTag("", ""))),
    posTrackToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("posTrackInputTag"))),
    negTrackToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("negTrackInputTag"))),
    genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleInputTag"))),
    outputCollectionName_(iConfig.getParameter<string>("outputCollectionName")),
    isSignal_(iConfig.getParameter<bool>("isSignal")),
    debug_(iConfig.getParameter<int>("debug")),
    cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
    dxymax_(cutSet_.getParameter<double>("dxymax")),
    dzmax_(cutSet_.getParameter<double>("dzmax")),
    dRmin_(cutSet_.getParameter<double>("dRmin")),
    dRmax_(cutSet_.getParameter<double>("dRmax")),
    trkPairdRmax_(cutSet_.getParameter<double>("trkpairdRmax")),
    phiPairMmin_(cutSet_.getParameter<double>("phipairMmin")),
    phiPairMmax_(cutSet_.getParameter<double>("phipairMmax")),
    bsPtMin_(cutSet_.getParameter<double>("ptmin"))
{
  produces<TkBsCandidateCollection>(outputCollectionName_);
  usesResource(TFileService::kSharedResource);
    cout << "--> L1BsMesonSelectionProducer:"
       << " dxymax: " << dxymax_
       << " dzmax: " << dzmax_
       << " dRmin: " << dRmin_
       << " dRmax: " << dRmax_
       << " trkPairdRmax: " << trkPairdRmax_
       << " phiPairMmin: " << phiPairMmin_
       << " phiPairMmax: " << phiPairMmax_
       << " bsPtMin: " << bsPtMin_
       << endl;

}

L1BsMesonSelectionProducer::~L1BsMesonSelectionProducer() {}
void L1BsMesonSelectionProducer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return;
  
  for (size_t i = 0; i < evCounters_.size(); ++i) evCounters_[i] = 0;
  evtCutFlow_   = fs->make<TH1D>("evtCutFlow", "",11,-0.5,10.5);
  ptAllPhiH_     = fs->make<TH1D>("ptPhi", "Phi Candidate pT", 200, 0.0, 30.0);
  etaAllPhiH_    = fs->make<TH1D>("etaPhi", "Phi Candidate eta", 200, -3.0, 3.0);
  phiAllPhiH_    = fs->make<TH1D>("phiPhi", "Phi Candidate phi", 200, -M_PI, M_PI);
  massAllPhiH_   = fs->make<TH1D>("massPhi", "Phi Candidate mass", 200, 0.99, 1.04);
  z0AllPhiH_     = fs->make<TH1D>("z0Phi", "Phi Candidate d0", 150, -15, 15);

  nPhiH_      = fs->make<TH1D>("nPhi", "Number of Phi meson candidates", 20, -0.5, 19.5);
  nBsH_      = fs->make<TH1D>("nBs", "Number of Bs  candidates", 20, -0.5, 19.5);
  dzPhiPairH_ = fs->make<TH1D>("dzPhiPair", "dz between all phi pairs", 100, -5.0, 5.0);
  drPhiPairH_ = fs->make<TH1D>("drPhiPair", "dR angle between all phi pairs", 200, 0.0, 2.5);
  bsMassH_    = fs->make<TH1D>("bsmass","Bs Mass", 200, 5.0, 5.8);
  bsMassFH_    = fs->make<TH1D>("bsmassF", "Bs Mass Final", 200, 5.0, 5.8);
  bsMassBinH_    = fs->make<TH1D>("bsmassBin","Bs Mass 410 bins", 6554, 5.0, 5.8);
  bsMassFBinH_    = fs->make<TH1D>("bsmassBinF", "Bs Mass Final 410 bins",6554, 5.0, 5.8);

  bsPtH_      = fs->make<TH1D>("bspt", "Bs pT", 200, 8.0, 28.0);
  bsEtaH_     = fs->make<TH1D>("bseta", "Bs eta", 200, -3.0, 3.0);
  bsPhiH_     = fs->make<TH1D>("bsphi", "Bs phi", 200, -M_PI, M_PI);

  drPhi1TrkPairH_ = fs->make<TH1D>("drPhi1TrkPair", "dr between the track pair (Phi1)", 100, 0, 0.2);
  drPhi2TrkPairH_ = fs->make<TH1D>("drPhi2TrkPair", "dr between the track pair (Phi2)", 100, 0, 0.2);
  
  trkPtH_[0] = fs->make<TH1D>("trk1Pt", "Highest pT track pT (Bs candidate)", 200, 0.0, 20.);
  trkPtH_[1] = fs->make<TH1D>("trk2Pt", "Second highest pT track pT (Bs candidate)", 200, 0.0, 10.);
  trkPtH_[2] = fs->make<TH1D>("trk3Pt", "Third highest pT track pT (Bs candidate)", 200, 0.0, 10.);
  trkPtH_[3] = fs->make<TH1D>("trk4Pt", "Lowest pT track pT (Bs candidate)", 200, 0.0, 10.);

  trkEtaH_[0] = fs->make<TH1D>("trk1Eta", "Highest pT track eta (Bs candidate)", 100, -3, 3);
  trkEtaH_[1] = fs->make<TH1D>("trk2Eta", "Second highest pT track eta (Bs candidate)", 100, -3, 3);
  trkEtaH_[2] = fs->make<TH1D>("trk3Eta", "Third highest pT track eta (Bs candidate)", 100, -3, 3);
  trkEtaH_[3] = fs->make<TH1D>("trk4Eta", "Lowest pT track eta (Bs candidate)", 100, -3, 3);

  trkPhiH_[0] = fs->make<TH1D>("trk1Phi", "Highest pT track phi (Bs candidate)", 100, -4, 4);
  trkPhiH_[1] = fs->make<TH1D>("trk2Phi", "Second highest pT track phi (Bs candidate)", 100, -4, 4);
  trkPhiH_[2] = fs->make<TH1D>("trk3Phi", "Third highest pT track phi (Bs candidate)", 100, -4, 4);
  trkPhiH_[3] = fs->make<TH1D>("trk4Phi", "Lowest pT track phi (Bs Candidate)", 100, -4, 4);
  
  trkZ0H_[0] = fs->make<TH1D>("trk1Z0", "Highest pT track z0 (Bs candidate)", 100, -20, 20);
  trkZ0H_[1] = fs->make<TH1D>("trk2Z0", "Second highest pT track z0 phi (Bs candidate)", 100, -20, 20);
  trkZ0H_[2] = fs->make<TH1D>("trk3Z0", "Third highest pT track z0 (Bs candidate)", 100, -20, 20);
  trkZ0H_[3] = fs->make<TH1D>("trk4Z0", "Lowest pT track z0 (Bs Candidate)", 100, -20, 20);

  trknStubH_[0] = fs->make<TH1D>("trk1nStub", "Highest pT track nStub (Bs candidate)", 11, -0.5, 10.5);
  trknStubH_[1] = fs->make<TH1D>("trk2nStub", "Second highest pT track nStub (Bs candidate)", 11, -0.5, 10.5);
  trknStubH_[2] = fs->make<TH1D>("trk3nStub", "Third highest pT track nStub (Bs candidate)", 11, -0.5, 10.5);
  trknStubH_[3] = fs->make<TH1D>("trk4nStub", "Lowest pT track nStub (Bs candidate)", 11, -0.5, 10.5);

  trknStubPSH_[0] = fs->make<TH1D>("trk1nStubPS", "Highest pT track PS nStub (Bs candidate)", 8, -0.5, 7.5);
  trknStubPSH_[1] = fs->make<TH1D>("trk2nStubPS", "Second highest pT track PS nStub (Bs candidate)", 8, -0.5, 7.5);
  trknStubPSH_[2] = fs->make<TH1D>("trk3nStubPS", "Third highest pT track PS nStub (Bs candidate)", 8, -0.5, 7.5);
  trknStubPSH_[3] = fs->make<TH1D>("trk4nStubPS", "Lowest pT track PS nStub (Bs candidate)", 8, -0.5, 7.5);

  trkChi2BendH_[0] = fs->make<TH1D>("trk1Chi2Bend", "Highest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);
  trkChi2BendH_[1] = fs->make<TH1D>("trk2Chi2Bend", "Second highest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);
  trkChi2BendH_[2] = fs->make<TH1D>("trk3Chi2Bend", "Third highest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);
  trkChi2BendH_[3] = fs->make<TH1D>("trk4Chi2Bend", "Lowest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);

  trkChi2ZH_[0] = fs->make<TH1D>("trk1Chi2Z", "Highest pT track reduced chi2Z (Bs candidate)", 100, 0, 10);
  trkChi2ZH_[1] = fs->make<TH1D>("trk2Chi2Z", "Second highest pT track reduced chi2Z (Bs candidate)", 100, 0, 10);
  trkChi2ZH_[2] = fs->make<TH1D>("trk3Chi2Z", "Third highest pT track reduced chi2Z (Bs candidate)", 100, 0, 10);
  trkChi2ZH_[3] = fs->make<TH1D>("trk4Chi2Z", "Lowest pT track reduced chi2Z (Bs candidate)", 100, 0, 10);

  trkChi2XYH_[0] = fs->make<TH1D>("trk1Chi2XY", "Highest pT track reduced chi2XY (Bs candidate)", 100, 0, 20);
  trkChi2XYH_[1] = fs->make<TH1D>("trk2Chi2XY", "Second highest pT track reduced chi2XY (Bs candidate)", 100, 0, 20);
  trkChi2XYH_[2] = fs->make<TH1D>("trk3Chi2XY", "Third highest pT track reduced chi2XY (Bs candidate)", 100, 0, 20);
  trkChi2XYH_[3] = fs->make<TH1D>("trk4Chi2XY", "Lowest pT track reduced chi2XY (Bs candidate)", 100, 0, 20);

  ptPhiH_[0]   = fs->make<TH1D>("ptPhi1",   "Phi1 pT (Bs Candidate)", 200, 4.0, 24.);
  etaPhiH_[0]  = fs->make<TH1D>("etaPhi1",  "Phi1 Eta (Bs Candidate)", 200, -3.0, 3.0);
  phiPhiH_[0]  = fs->make<TH1D>("phiPhi1",  "Phi1 Phi (Bs Candidate)", 200, -M_PI, M_PI);
  massPhiH_[0] = fs->make<TH1D>("massPhi1", "Phi1 mass (Bs Candidate)", 200, 0.98, 1.03);
  z0PhiH_[0]   = fs->make<TH1D>("z0Phi1",   "Phi1 Z0 (Bs candidate)", 200, -15.0, 15.0);

  ptPhiH_[1]   = fs->make<TH1D>("ptPhi2",   "Phi2 pT (Bs Candidate)", 200, 4.0, 24.);
  etaPhiH_[1]  = fs->make<TH1D>("etaPhi2",  "Phi2 Eta (Bs Candidate)", 200, -3.0, 3.0);
  phiPhiH_[1]  = fs->make<TH1D>("phiPhi2",  "Phi2 Phi (Bs Candidate)", 200, -M_PI, M_PI);
  massPhiH_[1] = fs->make<TH1D>("massPhi2", "Phi2 mass (Bs Candidate)", 200, 0.99, 1.04);
  z0PhiH_[1]   = fs->make<TH1D>("z0Phi2",   "Phi2 Z0 (Bs candidate)", 200, -15.0, 15.0);
  
  if (isSignal_) {
    drGenPhiH_   = fs->make<TH1D>("drGenPhi1H",   "dR GenPhi and Phi Candidate", 200, 0.0, 1.0);
    dptGenPhiH_  = fs->make<TH1D>("dptGenPhi1H",  "dPt GenPhi and Phi Candidate", 200, -5.0, 5.0);
    detaGenPhiH_ = fs->make<TH1D>("detaGenPhi1H", "dEta GenPhi and Phi Candidate", 40, 0.0, 5.0);
  }
}
// ------------ method called to produce the data  ------------
//void L1BsMesonSelectionProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
void L1BsMesonSelectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)  {

  auto L1BsMesonOutput = std::make_unique<TkBsCandidateCollection>();

  ++evCounters_[0];
  evtCutFlow_->Fill(0);

  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  bool res = false;
  if (isSignal_) res = iEvent.getByToken(genParticleToken_, genParticleHandle);

  TTTrackCollectionHandle posTrackHandle;
  iEvent.getByToken(posTrackToken_, posTrackHandle);

  TTTrackCollectionHandle negTrackHandle;
  iEvent.getByToken(negTrackToken_, negTrackHandle);

  if (posTrackHandle->size() < 2 || negTrackHandle->size() < 2) return;
  ++evCounters_[1];
  evtCutFlow_->Fill(1);

  // Tracker Topology
  const TrackerTopology& tTopo = iSetup.getData(tTopoToken_);

  TkPhiCandidateCollectionHandle l1PhiCandHandle;
  iEvent.getByToken(l1PhiCandToken_, l1PhiCandHandle);
  size_t nPhiMesonOutputApproximate = l1PhiCandHandle->size();

  for (size_t i = 0; i < nPhiMesonOutputApproximate; i++) {
    const auto& obj = l1PhiCandHandle->at(i);
    ptAllPhiH_->Fill(obj.p4().Pt());
    etaAllPhiH_->Fill(obj.p4().Eta());
    phiAllPhiH_->Fill(obj.p4().Phi());
    massAllPhiH_->Fill(obj.p4().M());
    z0AllPhiH_->Fill(obj.vz());
  }
  nPhiH_->Fill(nPhiMesonOutputApproximate);
  if(nPhiMesonOutputApproximate < 2) return;
  ++evCounters_[2];
  evtCutFlow_->Fill(2);
  size_t nBsMesonOutputApproximate = 2 * nPhiMesonOutputApproximate;
  L1BsMesonOutput->reserve(nBsMesonOutputApproximate);

  std::array<unsigned int, 8> iCounters = {{0, 0, 0, 0, 0, 0, 0,0}};

  
  for (size_t i = 0; i < nPhiMesonOutputApproximate; i++) {
    const auto& phiCand1 = l1PhiCandHandle->at(i);
    const math::XYZTLorentzVector& phiCand1P4 = phiCand1.p4();
    double ptPhi1   = phiCand1P4.Pt();
    double etaPhi1  = phiCand1P4.Eta();
    double phiPhi1  = phiCand1P4.Phi();
    double massPhi1 = phiCand1P4.M();       
    math::PtEtaPhiMLorentzVector phi1P4(ptPhi1, etaPhi1, phiPhi1, massPhi1);    

    const auto& trka_p = phiCand1.trkPtr(0);
    const auto& trka_n = phiCand1.trkPtr(1);
    double drTrkPhi1 = phiCand1.dRTrkPair();
      
    for (size_t j = i+1; j < nPhiMesonOutputApproximate; j++) {
      const auto& phiCand2 = l1PhiCandHandle->at(j);
      // Must ensure that the 2 Phi mesons are made of 4 distinct tracks      
      // this is non-trivial if the reconstructed Phi mesons do not store reference to the parent tracks
      const auto& trkb_p = phiCand2.trkPtr(0);
      const auto& trkb_n = phiCand2.trkPtr(1);
      
      if (trka_p == trkb_p || trka_n == trkb_n) continue;
      ++iCounters[0];

      if (duplicateTrk(trka_p, trkb_p, tTopo) || duplicateTrk(trka_n, trkb_n, tTopo)) continue;
      ++iCounters[1];
      
     
      const math::XYZTLorentzVector& phiCand2P4 = phiCand2.p4();
      double ptPhi2   = phiCand2P4.Pt();
      double etaPhi2  = phiCand2P4.Eta();
      double phiPhi2  = phiCand2P4.Phi();
      double massPhi2 = phiCand2P4.M();       
      math::PtEtaPhiMLorentzVector phi2P4(ptPhi2, etaPhi2, phiPhi2, massPhi2);    
      
      math::PtEtaPhiMLorentzVector tempP4 = phi1P4 + phi2P4;
      math::XYZTLorentzVector bsP4(tempP4.Px(), tempP4.Py(), tempP4.Pz(), tempP4.E());
      TkBsCandidate tkBs(bsP4, phiCand1, phiCand2);      
      //if (tkBs.dxyPhiPair() > dxymax_) continue;
      double dzPhiPair = tkBs.dzPhiPair();
      dzPhiPairH_->Fill(dzPhiPair);
      if (std::fabs(dzPhiPair) > dzmax_) continue;
      ++iCounters[2];
      
      double drPhiPair = tkBs.dRPhiPair();
      drPhiPairH_->Fill(drPhiPair);
      if (drPhiPair < dRmin_ || drPhiPair > dRmax_) continue;
      ++iCounters[3];

      double drTrkPhi2 = phiCand2.dRTrkPair();
      drPhi1TrkPairH_->Fill(drTrkPhi1);
      drPhi2TrkPairH_->Fill(drTrkPhi2);      

      if (drTrkPhi1 > trkPairdRmax_ || drTrkPhi2 > trkPairdRmax_) continue;
      ++iCounters[4];

      double mass = tkBs.p4().M();
      bsMassH_->Fill(mass);
      bsMassBinH_->Fill(mass);

      if (mass < phiPairMmin_ || mass > phiPairMmax_) continue;
      ++iCounters[5];

      bsPtH_->Fill(tkBs.p4().Pt());
      bsEtaH_->Fill(tkBs.p4().Eta());
      bsPhiH_->Fill(tkBs.p4().Phi());
      
      if (tkBs.p4().Pt() < bsPtMin_) continue;
      ++iCounters[6];

      //Gen Match
      if (isSignal_) {
        if (res && genParticleHandle.isValid()) {
          reco::GenParticleCollection genParticles = *genParticleHandle.product();
          if (!isGenMatched(genParticles, {phiCand1, phiCand2})) continue;
          ++iCounters[7];
        }
        else {
          cerr << "filter: GenParticleCollection for InputTag genParticles not found!" << endl;
        }
      }
      else 
	++iCounters[7];
      bsMassFH_->Fill(mass);
      bsMassFBinH_->Fill(mass);

      L1BsMesonOutput->push_back(tkBs);



      // Track properties
      vector<edm::Ptr<L1TTTrackType>> trackList {
        trka_p, trka_n, trkb_p, trkb_n
      };
      std::sort(std::begin(trackList), std::end(trackList), 
		[] (auto lhs, auto rhs) {
		  return lhs->momentum().perp() > rhs->momentum().perp();
		});

      for (size_t i = 0; i < trackList.size(); ++i) {
	fillTrackQuantities(trackList[i], tTopo,
			    trkPtH_[i], trkEtaH_[i], trkPhiH_[i], trkZ0H_[i], 
			    trknStubH_[i], trknStubPSH_[i], trkChi2BendH_[i], 
			    trkChi2ZH_[i], trkChi2XYH_[i]);
      }
      // Phi meson properties
      ptPhiH_[0]->Fill(ptPhi1);
      etaPhiH_[0]->Fill(etaPhi1);
      phiPhiH_[0]->Fill(phiPhi1);
      massPhiH_[0]->Fill(massPhi1);
      z0PhiH_[0]->Fill(phiCand1.vz());

      ptPhiH_[1]->Fill(ptPhi2);
      etaPhiH_[1]->Fill(etaPhi2);
      phiPhiH_[1]->Fill(phiPhi2);
      massPhiH_[1]->Fill(massPhi2);
      z0PhiH_[1]->Fill(phiCand2.vz());
    }
  }
  for (size_t i = 0; i < iCounters.size(); ++i) {
    if (iCounters[i] > 0) {
      ++evCounters_[3+i];
      evtCutFlow_->Fill(3+i);
    }
  }
  
  nBsH_->Fill(L1BsMesonOutput->size());
  // Put the outputs into the event
  iEvent.put(std::move(L1BsMesonOutput), outputCollectionName_);

}


void L1BsMesonSelectionProducer::fillTrackQuantities(edm::Ptr<L1TTTrackType> trk, 
						     const TrackerTopology& tTopo, 
						     TH1D* trkPtH, 
						     TH1D* trkEtaH, 
						     TH1D* trkPhiH, 
						     TH1D* trkZ0H, 
						     TH1D* trknStubH, 
						     TH1D* trknStubPSH,
						     TH1D* trkChi2BendH, 
						     TH1D* trkChi2ZH, 
						     TH1D* trkChi2XYH)
{
  trkPtH->Fill(trk->momentum().perp());
  trkEtaH->Fill(trk->momentum().eta());
  trkPhiH->Fill(trk->momentum().phi());
  trkZ0H->Fill(trk->z0());
  
  int nStub, nStubSS, nStubPS;
  stubInfo(trk, tTopo, nStub, nStubSS, nStubPS);
  trknStubH->Fill(nStub);
  trknStubPSH->Fill(nStubPS);
  
  trkChi2BendH->Fill(trk->stubPtConsistency());
  trkChi2ZH->Fill(trk->chi2ZRed());
  trkChi2XYH->Fill(trk->chi2XYRed());
}


bool L1BsMesonSelectionProducer::duplicateTrk(edm::Ptr<L1TTTrackType> trka,
					      edm::Ptr<L1TTTrackType> trkb,
					      const TrackerTopology& tTopo)
{
  int nStub_trka, nStubSS_trka, nStubPS_trka;
  stubInfo(trka, tTopo, nStub_trka, nStubSS_trka, nStubPS_trka);

  int nStub_trkb, nStubSS_trkb, nStubPS_trkb;
  stubInfo(trkb, tTopo, nStub_trkb, nStubSS_trkb, nStubPS_trkb);

  bool dupl = fabs(trka->momentum().eta()  - trkb->momentum().eta()) < 1.0e-03 &&
	      fabs(trka->momentum().phi()  - trkb->momentum().phi()) < 1.0e-03 &&
	      fabs(trka->momentum().perp() - trkb->momentum().perp()) < 3.0e-02 &&
              nStub_trka == nStub_trkb &&
              nStubPS_trka == nStubPS_trkb &&
	      fabs(trka->z0() - trkb->z0()) < 3.0e-01;
  if (dupl) {
#ifdef __DEBUG__
    printHeader(2);
    printTrk(trka, tTopo, 1);
    printTrk(trkb, tTopo, 2);
#endif
    return true;
  }
  return false;
}
void L1BsMesonSelectionProducer::printHeader(int ntrk) {
  cout << " Number of Tracks " << ntrk << endl
       << "  itrk       pT      Eta      Phi     chi2Bend  chi2RZ chi2RPhi"
       << " Curvature  VertexX  VertexY  VertexZ       d0       z0 nStub nStubPS"
       << endl
       << std::setiosflags(std::ios::fixed)
       << std::setprecision(2);
}
void L1BsMesonSelectionProducer::printTrk(edm::Ptr<L1TTTrackType> trkPtr, const TrackerTopology& tTopo, int itrk) {
  int nStub, nStubSS, nStubPS;
  stubInfo(trkPtr, tTopo, nStub, nStubSS, nStubPS);

  cout << setw(6) << itrk
       << std::setprecision(2)
       << setw(9) << trkPtr->momentum().perp()
       << setw(9) << trkPtr->momentum().eta()
       << setw(9) << trkPtr->momentum().phi()
       << setw(9) << trkPtr->stubPtConsistency()
       << std::setprecision(3)
       << setw(9) << trkPtr->chi2ZRed()
       << setw(9) << trkPtr->chi2XYRed()
       << std::setprecision(4)
       << setw(10) << trkPtr->rInv()
       << std::setprecision(3)
       << setw(9) << trkPtr->POCA().x()
       << setw(9) << trkPtr->POCA().y()
       << setw(9) << trkPtr->POCA().z()
       << setw(9) << trkPtr->d0()
       << setw(9) << trkPtr->z0()
       << setw(6) << nStub
       << setw(8) << nStubPS
       << endl;
}
void L1BsMesonSelectionProducer::stubInfo(edm::Ptr<L1TTTrackType> trk,
					  const TrackerTopology& tTopo,
					  int& nStub,
					  int& nStubSS,
					  int& nStubPS)
{
  L1TTStubCollection stubs = trk->getStubRefs();
  nStub = stubs.size();
  nStubPS = 0;
  nStubSS = 0;
  for (const auto& obj: stubs) {
    const DetId& detid = obj->getDetId();
    if (detid.det() != DetId::Detector::Tracker) continue;
    if (detid.subdetId() == StripSubdetector::TOB) {
      (tTopo.tobLayer(detid) <= 3) ? nStubPS++ : nStubSS++;
    }
    else if (detid.subdetId() == StripSubdetector::TID) {
      (tTopo.tidRing(detid) <= 9) ? nStubPS++ : nStubSS++;
    }
  }
}

void L1BsMesonSelectionProducer::endJob() {
  ostringstream bsPtCutTag;
  bsPtCutTag << ">= 1 Bs candidate with pT > " << fixed << setprecision(1) << bsPtMin_ << " GeV";

  ostringstream bsMassCutTag;
  bsMassCutTag << fixed << setprecision(1) << phiPairMmin_ << " < Bs Mass < " << phiPairMmax_ << " GeV";

  ostringstream dzPhiPairCutTag;
  dzPhiPairCutTag << "dzPhiPair < " << fixed << setprecision(1)
                  << dzmax_ << " cm";

  ostringstream drPhiPairCutTag;
  drPhiPairCutTag << fixed << setprecision(1) << dRmin_ << " < drPhiPair < " << dRmax_;

  ostringstream drPhiTrkCutTag;
  drPhiTrkCutTag << "dR(Phi Track pair) < " << fixed << setprecision(2) << trkPairdRmax_;

  std::vector<std::string> tags {
    "                Total",
    " ntrk_p>=2 && ntrk_n>=2",
      "             #Phi > 1",
      "        4 distinct tracks",
      "    duplicate track removal",
      dzPhiPairCutTag.str(),
      drPhiPairCutTag.str(),
      drPhiTrkCutTag.str(),
      bsMassCutTag.str(),
      bsPtCutTag.str(),
      "Gen matched",
      };

  /*  std::cout<<"---------------------------------------------"<<std::endl;
    std::cout<<"----------Cutflow Bs-Evt Simulation----------"<<std::endl;
    std::cout<<"---------------------------------------------"<<std::endl;
    unsigned int nlabels = (isSignal_) ? tags.size() : tags.size() -1;
    for (size_t i = 0; i < nlabels; ++i)
    std::cout << setw(40) << tags[i] << setw(10) << evCounters_[i] << std::endl;*/

  int nbins = evtCutFlow_->GetNbinsX();
  if (!isSignal_) --nbins;
  for (int i = 1; i <=nbins; ++i) {
    evtCutFlow_->GetXaxis()->SetBinLabel(i, tags[i-1].c_str());
  }
  cout << "==> L1BsMesonSelectionProducer::endJob()" << endl;
  unsigned int nlabels = (isSignal_) ? tags.size() : tags.size() -1;
  for (size_t i = 0; i < nlabels; ++i)
    cout << setw(40) << tags[i] << setw(10) << evCounters_[i] << endl;

}



bool L1BsMesonSelectionProducer::isGenMatched(const reco::GenParticleCollection& genParticles,
                                              const vector<TkPhiCandidate>& phiList,
                                              bool verbose) {
  vector<const reco::Candidate*> genPhiList;
  for (auto it = genParticles.begin(); it != genParticles.end(); ++it) {
    const reco::Candidate& genp = *it;
    if (abs(genp.pdgId()) != 333) continue; // Phi                                                                                                                                                                                                                                                                                                                            

    // Find mother                                                                                                                                                                                                                                                                                                                                                            
    reco::Candidate* m = const_cast<reco::Candidate*>(genp.mother());
    if (m == nullptr) continue;
    while (m->pdgId() == genp.pdgId()) {
      m = const_cast<reco::Candidate*>(m->mother());
    }
    if (m == nullptr || abs(m->pdgId()) != 531) continue; // Bs                                                                                                                                                                                                                                                                                                              \
                                                                                                                                                                                                                                                                                                                                                                              

    unsigned int ndau = 0;
    for (size_t j = 0; j < genp.numberOfDaughters(); ++j) {
      const reco::Candidate* d = genp.daughter(j);
      if (abs(d->pdgId()) != 321 || d->pt() < 2.0 || fabs(d->eta()) > 2.4) continue;
      ++ndau;
    }
    if (ndau < 2) continue;
    genPhiList.push_back(&genp);
  }
  if (genPhiList.size() != 2) return false;

  unsigned int nmatch = 0;
  for (const auto& recophi: phiList) {
    double ptRecoPhi   = recophi.pt();
    double etaRecoPhi  = recophi.eta();
    double phiRecoPhi  = recophi.phi();

    double min_dr{999.9};
    reco::Candidate* mGenPhi = nullptr;
    for (const reco::Candidate* genphi: genPhiList) {
      double ptGenPhi  = genphi->pt();
      double etaGenPhi = genphi->eta();
      double phiGenPhi = genphi->phi();

      double dr = sqrt(pow((phiRecoPhi - phiGenPhi), 2) + pow((etaRecoPhi - etaGenPhi), 2));
      if (dr < min_dr) {
        min_dr = dr;
        mGenPhi = const_cast<reco::Candidate*>(genphi);
      }
    }

    if (mGenPhi != nullptr) {
      drGenPhiH_->Fill(min_dr);
      dptGenPhiH_->Fill(ptRecoPhi - mGenPhi->pt());
      detaGenPhiH_->Fill(etaRecoPhi - mGenPhi->eta());
    }

    if (min_dr > 0.05) continue;
    ++nmatch;
  }
  if (nmatch == 2) return true;
  return false;
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1BsMesonSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PhiCandInputTag", edm::InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"));
  desc.add<edm::InputTag>("posTrackInputTag", edm::InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedPositivecharge"));
  desc.add<edm::InputTag>("negTrackInputTag", edm::InputTag("L1KaonTrackSelectionProducer","Level1TTKaonTracksSelectedNegativecharge"));
  desc.add<edm::InputTag>("GenParticleInputTag", edm::InputTag("genParticles", ""));
  desc.add<string>("outputCollectionName", "Level1TTKaonTracksSelected");
  desc.add<bool>("isSignal", true)->setComment("true for signal, false for Rate Estimate");
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dxymax", 999.0)->setComment("dxy must be less than this value, [cm]");
    descCutSet.add<double>("dzmax", 0.5)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("dRmin", 0.12)->setComment("dr must be greater than this value, []");
    descCutSet.add<double>("dRmax", 1.0)->setComment("dr must be less than this value, []");
    descCutSet.add<double>("trkpairdRmax", 0.2)->setComment("track pair opening angle must be less than this value, [GeV]");
    descCutSet.add<double>("phipairMmin", 5.0)->setComment("phipair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("phipairMmax", 5.8)->setComment("phipair mass must be less than this value, [GeV]");
    descCutSet.add<double>("ptmin", 13.0)->setComment("Bs pT must be greater than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);
  }
  descriptions.addWithDefaultLabel(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1BsMesonSelectionProducer);
