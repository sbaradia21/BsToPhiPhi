// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1BsMesonSelectionEmulationProducer
//
/**\class L1BsMesonSelectionProducer L1BsMesonSelectionProducer.cc L1Trigger/L1TTrackMatch/plugins/L1BsMesonSelectionProducer.cc                                  
 Description: Select Bs candidates from the reconstructed Phi candidates                                                                                          
 Implementation:                                                                                                                                       
     Input:                                                                                                                                             
         l1t::TkLightMesonWordCollection - A collection of reconstructed Phi candidates                                                                    
     Output:                                                                                                                                 
         l1t::TkLightMesonWordCollection - A collection of reconstructed Bs candidates                                                                        
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
#include <array>
// Xilinx HLS includes
#include <ap_fixed.h>
#include <ap_int.h>
#include <stdio.h>
#include <cassert>
#include <cstdlib>

// user include files
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidateFwd.h"
#include "DataFormats/L1Trigger/interface/TkLightMesonWord.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Gen-level stuff                                                                                                                                                                                                                
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "L1TrackWordUnpacker.h"
#include <TH1D.h>
#include <TH2D.h>

#include <TFile.h>

// class declaration
using namespace std;
using namespace edm;
using namespace l1t;

class L1BsMesonSelectionEmulationProducer : public edm::one::EDProducer<edm::one::SharedResources> {
//class L1BsMesonSelectionEmulationProducer : public edm::global::EDProducer<> {
public:
  using L1TTTrackType            = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollection        = std::vector<L1TTTrackType>;
  //using TTTrackRef               = edm::Ref<TTTrackCollection>; 
  using TTTrackRefCollection     = edm::RefVector<TTTrackCollection>;
  using TTTrackCollectionHandle  = edm::Handle<TTTrackRefCollection>;
  //using TTTrackRefCollectionUPtr = std::unique_ptr<TTTrackRefCollection>; 
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  explicit L1BsMesonSelectionEmulationProducer(const edm::ParameterSet&);
  ~L1BsMesonSelectionEmulationProducer();
  void beginJob() override;
  void endJob() override;
  static constexpr double KaonMass = 0.493677 ;
  static bool duplicateTrk(const L1TTTrackType& trka, const L1TTTrackType& trkb);

  static void fillTrackQuantities(const L1TTTrackType& trk,
                           TH1D* trkPtH,
                           TH1D* trkEtaH,
                           TH1D* trkPhiH,
                           TH1D* trkZ0H,
                           TH1D* trknStubH,
                           TH1D* trkChi2BendH,
                           TH1D* trkChi2RZH,
                           TH1D* trkChi2RPhiH);
  bool isGenMatched(const reco::GenParticleCollection& genParticles,
                    const vector<TkLightMesonWord>& phiList,
                    bool verbose=false);

  double ETAPHI_LSB = M_PI / (1 << 12);
  double Z0_LSB = 0.05;

private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word
  
  using TkLightMesonWordCollectionHandle = edm::Handle<TkLightMesonWordCollection>;

  //  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
  void produce(edm::Event&, const edm::EventSetup&) override;

  //void endJob();

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TkLightMesonWordCollection> l1PhiMesonWordToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> posTrackToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> negTrackToken_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  const std::string outputCollectionName_;
  
  const bool isSignal_;
  int debug_;
  const edm::ParameterSet cutSet_;
  const double dxymax_, dzmax_, dRmin_, dRmax_, trkPairdRmax_, phiPairMmin_, phiPairMmax_, bsPtMin_;

  std::array<unsigned int, 11> evCounters_;
  TH1D *evtCutFlow_;
  TH1D *ptPhiH_, *etaPhiH_, *phiPhiH_, *massPhiH_, *z0PhiH_;
  TH2D *ptpPhi2DH_, *ptnPhi2DH_;
  TH1D *nPhiH_, *nBsH_;
  TH1D *dzPhiPairH_, *drPhiPairH_, *drPhiPair2H_, *bsPtH_, *bsEtaH_, *bsPhiH_;
  TH1D *drPhi1TrackPairH_, *ptPhi1H_, *etaPhi1H_, *phiPhi1H_, *z0Phi1H_;
  TH1D *drPhi2TrackPairH_, *ptPhi2H_, *etaPhi2H_, *phiPhi2H_, *z0Phi2H_;

  TH1D *trk1PtH_,  *trk2PtH_,  *trk3PtH_,  *trk4PtH_;
  TH1D *trk1EtaH_, *trk2EtaH_, *trk3EtaH_, *trk4EtaH_;
  TH1D *trk1PhiH_, *trk2PhiH_, *trk3PhiH_, *trk4PhiH_;
  TH1D *trk1Z0H_,  *trk2Z0H_,  *trk3Z0H_,  *trk4Z0H_;
  TH1D *trk1nStubH_, *trk2nStubH_, *trk3nStubH_, *trk4nStubH_;
  TH1D *trk1Chi2BendH_, *trk2Chi2BendH_, *trk3Chi2BendH_, *trk4Chi2BendH_;
  TH1D *trk1Chi2RZH_, *trk2Chi2RZH_, *trk3Chi2RZH_, *trk4Chi2RZH_;
  TH1D *trk1Chi2RPhiH_, *trk2Chi2RPhiH_, *trk3Chi2RPhiH_, *trk4Chi2RPhiH_;

  TH1D *bsMassH_, *bsMass1H_, *bsMass1TrkH_, *bsMass1Trk1H_, *bsMassCheckH_;
  TH1D *phi1MassH_, *phi2MassH_;
  //  TH1D *bsMass2H_, *bsMass2TrkH_;
  TH1D *bsMassBinH_, *bsMass1BinH_;
  
  TH1D *drGenPhiH_, *dptGenPhiH_, *detaGenPhiH_;
  
};

//
// constructors and destructor
//
L1BsMesonSelectionEmulationProducer::L1BsMesonSelectionEmulationProducer(const edm::ParameterSet& iConfig)
  : l1PhiMesonWordToken_(consumes<TkLightMesonWordCollection>(iConfig.getParameter<edm::InputTag>("l1PhiMesonWordInputTag"))),
    posTrackToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("posTrackInputTag"))),
    negTrackToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("negTrackInputTag"))),
    //SB Gen Match
    genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleInputTag"))),
    outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
    isSignal_(iConfig.getParameter<bool>("isSignal")),
    debug_(iConfig.getParameter<int>("debug")),
    cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
    dxymax_(cutSet_.getParameter<double>("dxymax")),
    dzmax_(cutSet_.getParameter<double>("dzmax")),
    dRmin_(cutSet_.getParameter<double>("dRmin")),
    dRmax_(cutSet_.getParameter<double>("dRmax")),
    trkPairdRmax_(cutSet_.getParameter<double>("trkPairdRmax")),
    phiPairMmin_(cutSet_.getParameter<double>("phiPairMmin")),
    phiPairMmax_(cutSet_.getParameter<double>("phiPairMmax")),
    bsPtMin_(cutSet_.getParameter<double>("bsPtMin"))
{
  produces<l1t::TkLightMesonWordCollection>(outputCollectionName_);
}

L1BsMesonSelectionEmulationProducer::~L1BsMesonSelectionEmulationProducer() {}
void L1BsMesonSelectionEmulationProducer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return;
  for (size_t i = 0; i < evCounters_.size(); ++i) evCounters_[i] = 0;
  evtCutFlow_   = fs->make<TH1D>("evtCutFlow", "",11,-0.5,10.5);
  ptPhiH_     = fs->make<TH1D>("ptPhi", "Phi Candidate pT", 200, 0.0, 30.0);
  etaPhiH_    = fs->make<TH1D>("etaPhi", "Phi Candidate eta", 200, -3.0, 3.0);
  phiPhiH_    = fs->make<TH1D>("phiPhi", "Phi Candidate phi", 200, -M_PI, M_PI);
  massPhiH_   = fs->make<TH1D>("massPhi", "Phi Candidate mass", 200, 0.95, 1.15);
  z0PhiH_     = fs->make<TH1D>("z0Phi", "Phi Candidate d0", 200, -0.5, 0.5);

  ptpPhi2DH_  = fs->make<TH2D>("ptpPhi2D", "pTrk vs Phi Candidate pT", 200, 0.0, 100.0, 200, 0.0, 100.0);
  ptnPhi2DH_  = fs->make<TH2D>("ptnPhi2D", "nTrk vs Phi Candidate pT", 200, 0.0, 100.0, 200, 0.0, 100.0);

  
  nPhiH_      = fs->make<TH1D>("nPhi", "Number of Phi meson candidates", 20, -0.5, 19.5);
  nBsH_      = fs->make<TH1D>("nBs", "Number of Bs  candidates", 20, -0.5, 19.5);
  dzPhiPairH_   = fs->make<TH1D>("dzPhiPair", "dz(Phi1, Phi2)", 200, -0.5, 0.5);
  drPhiPairH_ = fs->make<TH1D>("drPhiPair", "dR angle between all phi pairs", 200, 0.0, 3.4);
  drPhiPair2H_ = fs->make<TH1D>("drPhiPair2", "dR angle between all phi pairs 2", 200, 0.0, 3.4);
  bsMassH_    = fs->make<TH1D>("bsmass", "Bs Mass after all cut from Bs word", 200, 5.0, 5.8);
  bsMassCheckH_    = fs->make<TH1D>("bsmassCheck", "Bs Mass after all cut from Bs word", 200, 5.0, 5.8);
  bsMassBinH_    = fs->make<TH1D>("bsmassBin", "Bs Mass after all cut from Bs word 410 bins",6554, 5.0, 5.8);

  bsMass1H_    = fs->make<TH1D>("bsmass1", "Bs Mass before mass cut", 200, 5.0, 5.8);
  bsMass1BinH_    = fs->make<TH1D>("bsmass1Bin", "Bs Mass before mass cut using Trks 410 bins", 6554, 5.0, 5.8);
  
  bsMass1TrkH_    = fs->make<TH1D>("bsmass1Trk", "Bs Mass before mass cut using Trks", 200, 5.0, 5.8);
  bsMass1Trk1H_   = fs->make<TH1D>("bsmass1Trk1", "Bs Mass before mass cut using Trks1", 200, 5.0, 5.8);

  //  bsMass1TrkBinH_    = fs->make<TH1D>("bsmass1TrkBin", "Bs Mass before mass cut using Trks 410 bins", 410, 5.0, 5.8);

  //  bsMass2H_    = fs->make<TH1D>("bsmass2", "Bs Mass after all cuts", 200, 5.0, 5.8);
  //  bsMass2TrkH_    = fs->make<TH1D>("bsmass2Trk", "Bs Mass after all cuts using Trks", 200, 5.0, 5.8);

  phi1MassH_     = fs->make<TH1D>("phi1Mass", "", 150, 1.0, 1.03);
  phi2MassH_     = fs->make<TH1D>("phi2Mass", "", 150, 1.0, 1.03);


  bsPtH_      = fs->make<TH1D>("bspt", "Bs pT", 200, 8.0, 28.0);
  bsEtaH_     = fs->make<TH1D>("bseta", "Bs eta", 200, -3.0, 3.0);
  bsPhiH_     = fs->make<TH1D>("bsphi", "Bs phi", 200, -M_PI, M_PI);


  drPhi1TrackPairH_ = fs->make<TH1D>("drPhi1TrackPair", "dR between the track pair (Phi1)", 100, 0, 0.2);
  drPhi2TrackPairH_ = fs->make<TH1D>("drPhi2TrackPair", "dR between the track pair (Phi2)", 100, 0, 0.2);

  ptPhi1H_  = fs->make<TH1D>("ptPhi1", "Phi 1 Candidate pT", 200, 0.0, 30.0);
  etaPhi1H_ = fs->make<TH1D>("etaPhi1", "Phi 1 Candidate eta", 200, -3.0, 3.0);
  phiPhi1H_ = fs->make<TH1D>("phiPhi1", "Phi 1 Candidate phi", 200, -M_PI, M_PI);
  z0Phi1H_  = fs->make<TH1D>("z0Phi1", "Phi 1 Candidate z0", 200, -15.0, 15.0);

  ptPhi2H_  = fs->make<TH1D>("ptPhi2", "Phi 2 Candidate pT", 200, 0.0, 30.0);
  etaPhi2H_ = fs->make<TH1D>("etaPhi2", "Phi 2 Candidate eta", 200, -3.0, 3.0);
  phiPhi2H_ = fs->make<TH1D>("phiPhi2", "Phi 2 Candidate phi", 200, -M_PI, M_PI);
  z0Phi2H_  = fs->make<TH1D>("z0Phi2", "Phi 2 Candidate z0", 200, -15.0, 15.0);

  trk1PtH_ = fs->make<TH1D>("trk1Pt", "Highest pT track pT (Bs candidate)", 200, 0.0, 20.);
  trk2PtH_ = fs->make<TH1D>("trk2Pt", "Second highest pT track pT (Bs candidate)", 200, 0.0, 10.);
  trk3PtH_ = fs->make<TH1D>("trk3Pt", "Third highest pT track pT (Bs candidate)", 200, 0.0, 10.);
  trk4PtH_ = fs->make<TH1D>("trk4Pt", "Lowest pT track pT (Bs candidate)", 200, 0.0, 10.);

  trk1EtaH_ = fs->make<TH1D>("trk1Eta", "Highest pT track eta (Bs candidate)", 100, -3, 3);
  trk2EtaH_ = fs->make<TH1D>("trk2Eta", "Second highest pT track eta (Bs candidate)", 100, -3, 3);
  trk3EtaH_ = fs->make<TH1D>("trk3Eta", "Third highest pT track eta (Bs candidate)", 100, -3, 3);
  trk4EtaH_ = fs->make<TH1D>("trk4Eta", "Lowest pT track eta (Bs candidate)", 100, -3, 3);

  trk1PhiH_ = fs->make<TH1D>("trk1Phi", "Highest pT track phi (Bs candidate)", 100, -4, 4);
  trk2PhiH_ = fs->make<TH1D>("trk2Phi", "Second highest pT track phi (Bs candidate)", 100, -4, 4);
  trk3PhiH_ = fs->make<TH1D>("trk3Phi", "Third highest pT track phi (Bs candidate)", 100, -4, 4);
  trk4PhiH_ = fs->make<TH1D>("trk4Phi", "Lowest pT track phi (Bs Candidate)", 100, -4, 4);

  trk1Z0H_ = fs->make<TH1D>("trk1Z0", "Highest pT track z0 (Bs candidate)", 100, -20, 20);
  trk2Z0H_ = fs->make<TH1D>("trk2Z0", "Second highest pT track z0 phi (Bs candidate)", 100, -20, 20);
  trk3Z0H_ = fs->make<TH1D>("trk3Z0", "Third highest pT track z0 (Bs candidate)", 100, -20, 20);
  trk4Z0H_ = fs->make<TH1D>("trk4Z0", "Lowest pT track z0 (Bs Candidate)", 100, -20, 20);

  trk1nStubH_ = fs->make<TH1D>("trk1nStub", "Highest pT track nStub (Bs candidate)", 11, -0.5, 10.5);
  trk2nStubH_ = fs->make<TH1D>("trk2nStub", "Second highest pT track nStub (Bs candidate)", 11, -0.5, 10.5);
  trk3nStubH_ = fs->make<TH1D>("trk3nStub", "Third highest pT track nStub (Bs candidate)", 11, -0.5, 10.5);
  trk4nStubH_ = fs->make<TH1D>("trk4nStub", "Lowest pT track nStub (Bs candidate)", 11, -0.5, 10.5);

  trk1Chi2BendH_ = fs->make<TH1D>("trk1Chi2Bend", "Highest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);
  trk2Chi2BendH_ = fs->make<TH1D>("trk2Chi2Bend", "Second highest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);
  trk3Chi2BendH_ = fs->make<TH1D>("trk3Chi2Bend", "Third highest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);
  trk4Chi2BendH_ = fs->make<TH1D>("trk4Chi2Bend", "Lowest pT track reduced chi2Bend (Bs candidate)", 100, 0, 5);

  trk1Chi2RZH_ = fs->make<TH1D>("trk1Chi2RZ", "Highest pT track reduced chi2RZ (Bs candidate)", 100, 0, 10);
  trk2Chi2RZH_ = fs->make<TH1D>("trk2Chi2RZ", "Second highest pT track reduced chi2RZ (Bs candidate)", 100, 0, 10);
  trk3Chi2RZH_ = fs->make<TH1D>("trk3Chi2RZ", "Third highest pT track reduced chi2RZ (Bs candidate)", 100, 0, 10);
  trk4Chi2RZH_ = fs->make<TH1D>("trk4Chi2RZ", "Lowest pT track reduced chi2RZ (Bs candidate)", 100, 0, 10);

  trk1Chi2RPhiH_ = fs->make<TH1D>("trk1Chi2RPhi", "Highest pT track reduced chi2RPhi (Bs candidate)", 100, 0, 20);
  trk2Chi2RPhiH_ = fs->make<TH1D>("trk2Chi2RPhi", "Second highest pT track reduced chi2RPhi (Bs candidate)", 100, 0, 20);
  trk3Chi2RPhiH_ = fs->make<TH1D>("trk3Chi2RPhi", "Third highest pT track reduced chi2RPhi (Bs candidate)", 100, 0, 20);
  trk4Chi2RPhiH_ = fs->make<TH1D>("trk4Chi2RPhi", "Lowest pT track reduced chi2RPhi (Bs candidate)", 100, 0, 20);
  if (isSignal_) {
    drGenPhiH_   = fs->make<TH1D>("drGenPhi1H",   "dR GenPhi and Phi Candidate", 200, 0.0, 1.0);
    dptGenPhiH_  = fs->make<TH1D>("dptGenPhi1H",  "dPt GenPhi and Phi Candidate", 200, -5.0, 5.0);
    detaGenPhiH_ = fs->make<TH1D>("detaGenPhi1H", "dEta GenPhi and Phi Candidate", 40, 0.0, 5.0);
  }
}

//
// member functions
//
// ------------ method called to produce the data  ------------
//void L1BsMesonSelectionEmulationProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
void L1BsMesonSelectionEmulationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)  {
  auto L1BsMesonEmulationOutput = std::make_unique<l1t::TkLightMesonWordCollection>();
  ++evCounters_[0];
  evtCutFlow_->Fill(0);

  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  bool res = false;
  if (isSignal_) res = iEvent.getByToken(genParticleToken_, genParticleHandle);

  TTTrackCollectionHandle posTrackHandle;
  iEvent.getByToken(posTrackToken_, posTrackHandle);
  
  TTTrackCollectionHandle negTrackHandle;
  iEvent.getByToken(negTrackToken_, negTrackHandle);

  if ((posTrackHandle->size() < 2) || (negTrackHandle->size() < 2)) return;
  ++evCounters_[1];
  evtCutFlow_->Fill(1);

  TkLightMesonWordCollectionHandle l1PhiMesonWordHandle;
  iEvent.getByToken(l1PhiMesonWordToken_, l1PhiMesonWordHandle);

  size_t nPhiMesonOutputApproximate = l1PhiMesonWordHandle->size();
  for (size_t i = 0; i < nPhiMesonOutputApproximate; i++) {
    const auto& obj = l1PhiMesonWordHandle->at(i);
    ptPhiH_->Fill(obj.pt());
    etaPhiH_->Fill(obj.glbeta());
    phiPhiH_->Fill(obj.glbphi());
    massPhiH_->Fill(obj.mass());
    z0PhiH_->Fill(obj.z0());

    unsigned int pIndexPhi = obj.firstIndex();
    unsigned int nIndexPhi = obj.secondIndex();
    const auto& pTrackPhi = *(posTrackHandle->at(pIndexPhi));
    const auto& nTrackPhi = *(negTrackHandle->at(nIndexPhi));
    double pt_pTrackPhi = l1trackunpacker::FloatPtFromBits(pTrackPhi), pt_nTrackPhi = l1trackunpacker::FloatPtFromBits(nTrackPhi);
    ptpPhi2DH_->Fill(pt_pTrackPhi,obj.pt());
    ptnPhi2DH_->Fill(pt_nTrackPhi,obj.pt());

  }
  nPhiH_->Fill(nPhiMesonOutputApproximate);
  if(nPhiMesonOutputApproximate < 2) return;
  ++evCounters_[2];
  evtCutFlow_->Fill(2);

  size_t nBsMesonOutputApproximate = 2 * nPhiMesonOutputApproximate;

  L1BsMesonEmulationOutput->reserve(nBsMesonOutputApproximate);

  std::array<unsigned int, 8> iCounters = {{0, 0, 0, 0, 0, 0, 0, 0}};
  for (size_t i = 0; i < nPhiMesonOutputApproximate; i++) {
    const auto& phiMesonWord1 = l1PhiMesonWordHandle->at(i);        
    double ptPhi1   = phiMesonWord1.pt();
    double etaPhi1  = phiMesonWord1.glbeta();
    double phiPhi1  = phiMesonWord1.glbphi();
    double z0Phi1   = phiMesonWord1.z0();
    double massPhi1 = 1.019;//phiMesonWord1.mass();
    unsigned int pIndexPhi1 = phiMesonWord1.firstIndex();
    unsigned int nIndexPhi1 = phiMesonWord1.secondIndex();
    const auto& pTrackPhi1 = *(posTrackHandle->at(pIndexPhi1));
    const auto& nTrackPhi1 = *(negTrackHandle->at(nIndexPhi1));

    double pxPhi1 = ptPhi1 * cos(phiPhi1);
    double pyPhi1 = ptPhi1 * sin(phiPhi1);
    double pzPhi1 = ptPhi1 * sinh(etaPhi1);
    double ePhi1  = sqrt(pow(ptPhi1, 2) + pow(pzPhi1, 2) + pow(massPhi1, 2));    
    math::PtEtaPhiMLorentzVector phi1P4(ptPhi1, etaPhi1, phiPhi1, massPhi1);    

    for (size_t j = i+1; j < nPhiMesonOutputApproximate; j++) {
      const auto& phiMesonWord2 = l1PhiMesonWordHandle->at(j);
      
      double ptPhi2   = phiMesonWord2.pt();
      double etaPhi2  = phiMesonWord2.glbeta();
      double phiPhi2  = phiMesonWord2.glbphi();
      double z0Phi2   = phiMesonWord2.z0();
      double massPhi2 = 1.019;//phiMesonWord2.mass();

      // Ensure that the 2 Phi mesons are made of 4 distinct tracks
      unsigned int pIndexPhi2 = phiMesonWord2.firstIndex();
      unsigned int nIndexPhi2 = phiMesonWord2.secondIndex();
      if (pIndexPhi1 == pIndexPhi2 || nIndexPhi1 == nIndexPhi2) continue;
      ++iCounters[0];

      const auto& pTrackPhi2 = *(posTrackHandle->at(pIndexPhi2));
      const auto& nTrackPhi2 = *(negTrackHandle->at(nIndexPhi2));
      if (duplicateTrk(pTrackPhi1, pTrackPhi2) || duplicateTrk(nTrackPhi1, nTrackPhi2)) continue;
      ++iCounters[1];

      double pxPhi2 = ptPhi2 * cos(phiPhi2);
      double pyPhi2 = ptPhi2 * sin(phiPhi2);
      double pzPhi2 = ptPhi2 * sinh(etaPhi2);
      double ePhi2  = sqrt(pow(ptPhi2, 2) + pow(pzPhi2, 2) + pow(massPhi2, 2));
      math::PtEtaPhiMLorentzVector phi2P4(ptPhi2, etaPhi2, phiPhi2, massPhi2);

      dzPhiPairH_->Fill(z0Phi1 - z0Phi2);
      if (std::fabs(z0Phi1 - z0Phi2) > dzmax_) continue;
      ++iCounters[2];

      double drPhiPair = reco::deltaR(etaPhi1, phiPhi1, etaPhi2, phiPhi2);
      double drPhiPair2 = sqrt(pow((phiPhi1 - phiPhi2), 2) + pow((etaPhi1 - etaPhi2), 2));
      drPhiPairH_->Fill(drPhiPair);
      drPhiPair2H_->Fill(drPhiPair2);
      if (drPhiPair < dRmin_ || drPhiPair > dRmax_) continue; 
      ++iCounters[3];

      // drTrackPair                                                                                                                                       
      double drTrackPairPhi1 = reco::deltaR(l1trackunpacker::FloatEtaFromBits(pTrackPhi1),
                                            l1trackunpacker::FloatPhiFromBits(pTrackPhi1),
                                            l1trackunpacker::FloatEtaFromBits(nTrackPhi1),
                                            l1trackunpacker::FloatPhiFromBits(nTrackPhi1));
      double drTrackPairPhi2 = reco::deltaR(l1trackunpacker::FloatEtaFromBits(pTrackPhi2),
                                            l1trackunpacker::FloatPhiFromBits(pTrackPhi2),
                                            l1trackunpacker::FloatEtaFromBits(nTrackPhi2),
                                            l1trackunpacker::FloatPhiFromBits(nTrackPhi2));
      drPhi1TrackPairH_->Fill(drTrackPairPhi1);
      drPhi2TrackPairH_->Fill(drTrackPairPhi2);

      if (drTrackPairPhi1 > trkPairdRmax_ ||  drTrackPairPhi2 > trkPairdRmax_) continue;
      ++iCounters[4];


      //double massPhiPair = sqrt(2 * ptPhi1 * ptPhi2 * (cosh(etaPhi1 - etaPhi2) - cos(phiPhi1 - phiPhi2)));
      double pxBs = pxPhi1 + pxPhi2;
      double pyBs = pyPhi1 + pyPhi2;
      double pzBs = pzPhi1 + pzPhi2;
      double eBs  = ePhi1  + ePhi2;

      double trk1pt=l1trackunpacker::FloatPtFromBits(pTrackPhi1);
      double trk2pt=l1trackunpacker::FloatPtFromBits(pTrackPhi2);
      double trk3pt=l1trackunpacker::FloatPtFromBits(nTrackPhi1);
      double trk4pt=l1trackunpacker::FloatPtFromBits(nTrackPhi2);

      double trk1eta=l1trackunpacker::FloatEtaFromBits(pTrackPhi1);
      double trk2eta=l1trackunpacker::FloatEtaFromBits(pTrackPhi2);
      double trk3eta=l1trackunpacker::FloatEtaFromBits(nTrackPhi1);
      double trk4eta=l1trackunpacker::FloatEtaFromBits(nTrackPhi2);

      double trk1phi=l1trackunpacker::FloatPhiFromBits(pTrackPhi1);
      double trk2phi=l1trackunpacker::FloatPhiFromBits(pTrackPhi2);
      double trk3phi=l1trackunpacker::FloatPhiFromBits(nTrackPhi1);
      double trk4phi=l1trackunpacker::FloatPhiFromBits(nTrackPhi2);

      double trk1px = trk1pt * cos(trk1phi), trk2px = trk2pt * cos(trk2phi), trk3px = trk3pt * cos(trk3phi), trk4px = trk4pt * cos(trk4phi);
      double trk1py = trk1pt * sin(trk1phi), trk2py = trk2pt * sin(trk2phi), trk3py = trk3pt * sin(trk3phi), trk4py = trk4pt * sin(trk4phi);
      double trk1pz = trk1pt * sinh(trk1eta), trk2pz = trk2pt * sinh(trk2eta), trk3pz = trk3pt * sinh(trk3eta), trk4pz = trk4pt * sinh(trk4eta);
      double trk1e = sqrt(pow(trk1pt,2) + pow(trk1pz,2) + pow(KaonMass,2)), trk2e = sqrt(pow(trk2pt,2) + pow(trk2pz,2) + pow(KaonMass,2)),trk3e = sqrt(pow(trk3pt,2) + pow(trk3pz,2) + pow(KaonMass,2)),trk4e = sqrt(pow(trk4pt,2) + pow(trk4pz,2) + pow(KaonMass,2));

      double pxBs_ = trk1px + trk2px + trk3px + trk4px;
      double pyBs_ = trk1py + trk2py + trk3py +	trk4py;
      double pzBs_ = trk1pz + trk2pz + trk3pz + trk4pz;
      double eBs_  = trk1e + trk2e + trk3e + trk4e;

      
      math::PtEtaPhiMLorentzVectorD kaon1_vec(trk1pt, trk1eta, trk1phi, KaonMass);
      math::PtEtaPhiMLorentzVectorD kaon2_vec(trk2pt, trk2eta, trk2phi, KaonMass);
      math::PtEtaPhiMLorentzVectorD kaon3_vec(trk3pt, trk3eta, trk3phi, KaonMass);
      math::PtEtaPhiMLorentzVectorD kaon4_vec(trk4pt, trk4eta, trk4phi, KaonMass);
      math::PtEtaPhiMLorentzVectorD phi1_vec = kaon1_vec+kaon3_vec;
      math::PtEtaPhiMLorentzVectorD phi2_vec = kaon2_vec + kaon4_vec;

      double massPhiPair = sqrt(pow(eBs, 2) - pow(pxBs, 2) - pow(pyBs, 2) - pow(pzBs, 2));
      //      double massPhiPairTrk = (phi1_vec+phi2_vec).M();
      double massPhiPairTrk = sqrt(pow(eBs_, 2) - pow(pxBs_, 2) - pow(pyBs_, 2) - pow(pzBs_, 2));
      double massPhiPairTrk1 = (kaon1_vec+kaon2_vec+kaon3_vec+kaon4_vec).M();
      
      //      double massPhiPairTrk = sqrt(pow(eBs_, 2) - pow(pxBs_, 2) - pow(pyBs_, 2) - pow(pzBs_, 2));
      
      bsMass1H_->Fill(massPhiPair);
      bsMass1TrkH_->Fill(massPhiPairTrk);
      bsMass1Trk1H_->Fill(massPhiPairTrk1);

      phi1MassH_ ->Fill(phi1_vec.M());
      phi2MassH_ ->Fill(phi2_vec.M());

      bsMass1BinH_->Fill(massPhiPair);
      //bsMass1TrkBinH_->Fill(massPhiPairTrk);
      
      
      //      bsMass1H_->Fill((phi1P4+phi2P4).M());
      
      if (massPhiPair < phiPairMmin_ || massPhiPair > phiPairMmax_) continue;
      ++iCounters[5];
      
      double pt = sqrt(pow(pxBs, 2) + pow(pyBs, 2));
      
      
      bsPtH_->Fill(pt);
      bsEtaH_->Fill(asinh(pzBs/sqrt(pow(pxBs, 2) + pow(pyBs, 2))));
      bsPhiH_->Fill(atan2(pyBs, pxBs));
      
      if (pt < bsPtMin_) continue;
      ++iCounters[6];

      if (isSignal_) {
        if (res && genParticleHandle.isValid()) {
          reco::GenParticleCollection genParticles = *genParticleHandle.product();
          if (!isGenMatched(genParticles, {phiMesonWord1, phiMesonWord2})) continue;
          ++iCounters[7];
        }
        else {
          cerr << "filter: GenParticleCollection for InputTag genParticles not found!" << endl;
        }
      }
      else
        ++iCounters[7];

      ptPhi1H_->Fill(ptPhi1);
      etaPhi1H_->Fill(etaPhi1);
      phiPhi1H_->Fill(phiPhi1);
      z0Phi1H_->Fill(z0Phi1);

      ptPhi2H_->Fill(ptPhi2);
      etaPhi2H_->Fill(etaPhi2);
      phiPhi2H_->Fill(phiPhi2);
      z0Phi2H_->Fill(z0Phi2);

      vector<L1TTTrackType> trackList {pTrackPhi1, nTrackPhi1, pTrackPhi2, nTrackPhi2};
      std::sort(std::begin(trackList), std::end(trackList),
                [] (const auto& lhs, const auto& rhs) {
                  return l1trackunpacker::FloatPtFromBits(lhs) > l1trackunpacker::FloatPtFromBits(rhs);
                });

      fillTrackQuantities(trackList[0], trk1PtH_, trk1EtaH_, trk1PhiH_, trk1Z0H_, trk1nStubH_, trk1Chi2BendH_, trk1Chi2RZH_, trk1Chi2RPhiH_);
      fillTrackQuantities(trackList[1], trk2PtH_, trk2EtaH_, trk2PhiH_, trk2Z0H_, trk2nStubH_, trk2Chi2BendH_, trk2Chi2RZH_, trk2Chi2RPhiH_);
      fillTrackQuantities(trackList[2], trk3PtH_, trk3EtaH_, trk3PhiH_, trk3Z0H_, trk3nStubH_, trk3Chi2BendH_, trk3Chi2RZH_, trk3Chi2RPhiH_);
      fillTrackQuantities(trackList[3], trk4PtH_, trk4EtaH_, trk4PhiH_, trk4Z0H_, trk4nStubH_, trk4Chi2BendH_, trk4Chi2RZH_, trk4Chi2RPhiH_);

      l1t::TkLightMesonWord::valid_t validBs = phiMesonWord1.valid() && phiMesonWord2.valid();
      l1t::TkLightMesonWord::pt_t ptBs = sqrt(pow(pxBs, 2) + pow(pyBs, 2));
      l1t::TkLightMesonWord::glbphi_t phiBs = atan2(pyBs, pxBs) / ETAPHI_LSB;
      l1t::TkLightMesonWord::glbeta_t etaBs = asinh(pzBs/sqrt(pow(pxBs, 2) + pow(pyBs, 2))) / ETAPHI_LSB;
      l1t::TkLightMesonWord::z0_t z0Bs = ((z0Phi1 + z0Phi2) / Z0_LSB) * 0.5;
      l1t::TkLightMesonWord::mass_t massBs  = sqrt(pow(eBs, 2) - pow(pxBs, 2) - pow(pyBs, 2) - pow(pzBs, 2));
      l1t::TkLightMesonWord::type_t typeBs = l1t::TkLightMesonWord::TkLightMesonTypes::kBsType;
      l1t::TkLightMesonWord::ntracks_t ntracksBs = 3;
      l1t::TkLightMesonWord::index_t firstPhiIndex     = i;
      l1t::TkLightMesonWord::index_t secondPhiIndex    = j;
      l1t::TkLightMesonWord::unassigned_t unassignedBs = 0;

      l1t::TkLightMesonWord bsWord(validBs, ptBs, phiBs, etaBs, z0Bs, massBs, typeBs, ntracksBs, firstPhiIndex, secondPhiIndex, unassignedBs); 
      //      l1t::TkLightMesonWord bsWord(validBs, ptBs, phiBs, etaBs, z0Bs, massBs, typeBs, ntracksBs, unassignedBs);      
      //      bsMassH_->Fill((phi1P4+phi2P4).M());
      bsMassH_->Fill(massBs.to_double());
      bsMassBinH_->Fill(massBs.to_double());
      
      //      bsMass2H_->Fill(massPhiPair);
      //bsMass2TrkH_->Fill(massPhiPairTrk);

      L1BsMesonEmulationOutput->push_back(bsWord);
    }
  }
  for (size_t i = 0; i < iCounters.size(); ++i) {
    if (iCounters[i] > 0) {
      ++evCounters_[3+i];
      evtCutFlow_->Fill(3+i);
    }
  }
  nBsH_->Fill(L1BsMesonEmulationOutput->size());
  if(L1BsMesonEmulationOutput->size()>0) {
    const auto& obj = L1BsMesonEmulationOutput->at(0);
    bsMassCheckH_->Fill(obj.mass());
  }
  // Put the outputs into the event
  iEvent.put(std::move(L1BsMesonEmulationOutput), outputCollectionName_);
}

bool L1BsMesonSelectionEmulationProducer::duplicateTrk(const L1TTTrackType& trka,
                                                       const L1TTTrackType& trkb)
{
  double eta_trka = l1trackunpacker::FloatEtaFromBits(trka),
    eta_trkb = l1trackunpacker::FloatEtaFromBits(trkb);
  double phi_trka = l1trackunpacker::FloatPhiFromBits(trka),
    phi_trkb = l1trackunpacker::FloatPhiFromBits(trkb);
  double pt_trka = l1trackunpacker::FloatPtFromBits(trka),
    pt_trkb = l1trackunpacker::FloatPtFromBits(trkb);
  double z0_trka = l1trackunpacker::FloatZ0FromBits(trka),
    z0_trkb = l1trackunpacker::FloatZ0FromBits(trkb);
  int nStub_trka = trka.getNStubs(),
    nStub_trkb = trkb.getNStubs();

  bool dupl = fabs(eta_trka - eta_trkb) < 1.0e-03 &&
	      fabs(phi_trka - phi_trkb) < 1.0e-03 &&
	      fabs(pt_trka  - pt_trkb)  < 3.0e-02 &&
              nStub_trka == nStub_trkb &&
	      fabs(z0_trka - z0_trkb) < 3.0e-01;
  if (dupl) return true;
  return false;
}
void L1BsMesonSelectionEmulationProducer::fillTrackQuantities(const L1TTTrackType& trk,
                                                              TH1D* trkPtH,
                                                              TH1D* trkEtaH,
                                                              TH1D* trkPhiH,
                                                              TH1D* trkZ0H,
                                                              TH1D* trknStubH,
                                                              TH1D* trkChi2BendH,
                                                              TH1D* trkChi2RZH,
                                                              TH1D* trkChi2RPhiH)
{
  trkEtaH->Fill(l1trackunpacker::FloatEtaFromBits(trk));
  trkPhiH->Fill(l1trackunpacker::FloatPhiFromBits(trk));
  trkPtH->Fill(l1trackunpacker::FloatPtFromBits(trk));
  trkZ0H->Fill(l1trackunpacker::FloatZ0FromBits(trk));
  trknStubH->Fill(trk.getNStubs());
  trkChi2BendH->Fill(trk.getBendChi2());
  trkChi2RZH->Fill(trk.getChi2RZ());
  trkChi2RPhiH->Fill(trk.getChi2RPhi());
}

void L1BsMesonSelectionEmulationProducer::endJob() {
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
      "            ntrk_p>=2 && ntrk_n>=2",
      "             #Phi > 1",
      "        4 distinct tracks",
      " duplicate track removal",
      dzPhiPairCutTag.str(),
      drPhiPairCutTag.str(),
      drPhiTrkCutTag.str(),
      bsMassCutTag.str(),
      bsPtCutTag.str(),
      "Gen  Matched"
      };
  int nbins = evtCutFlow_->GetNbinsX();
  if (!isSignal_) --nbins;
  for (int i = 1; i <=nbins; ++i) {
    evtCutFlow_->GetXaxis()->SetBinLabel(i, tags[i-1].c_str());
  }
  cout << "==> L1BsMesonSelectionEmulationProducer::endJob()" << endl;
  unsigned int nlabels = (isSignal_) ? tags.size() : tags.size() -1;
  for (size_t i = 0; i < nlabels; ++i)
    cout << setw(40) << tags[i] << setw(10) << evCounters_[i] << endl;
}

bool L1BsMesonSelectionEmulationProducer::isGenMatched(const reco::GenParticleCollection& genParticles,
                                                       const vector<TkLightMesonWord>& phiList,
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
    if (m == nullptr || abs(m->pdgId()) != 531) continue; // Bs                                                                                                                                                                                                                                                                                                                                                                                                           

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
    double etaRecoPhi  = recophi.glbeta();
    double phiRecoPhi  = recophi.glbphi();

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
void L1BsMesonSelectionEmulationProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PhiMesonWordInputTag", edm::InputTag("L1PhiMesonSelectionEmulationProducer","Level1TTPhiMesonSelectedEmulation"));
  desc.add<edm::InputTag>("posTrackInputTag", edm::InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedEmulationPositivecharge"));
  desc.add<edm::InputTag>("negTrackInputTag", edm::InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedEmulationNegativecharge"));
  desc.add<edm::InputTag>("GenParticleInputTag", edm::InputTag("genParticles", ""));
  desc.add<std::string>("outputCollectionName", "Level1TTBsMesonSelectedEmulation");
  desc.add<bool>("isSignal", true)->setComment("true for signal, false for Rate Estimate");
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dxymax", 999.0)->setComment("dxy must be less than this value, [cm]");
    descCutSet.add<double>("dzmax", 0.5)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("dRmin", 0.12)->setComment("dr must be greater than this value, []");
    descCutSet.add<double>("dRmax", 1.0)->setComment("dr must be less than this value, []");
    descCutSet.add<double>("trkPairdRmax", 0.2)->setComment("dR must be less than this value, [cm]");
    descCutSet.add<double>("phiPairMmin", 5.0)->setComment("phiPair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("phiPairMmax", 5.8)->setComment("phiPair mass must be less than this value, [GeV]");
    descCutSet.add<double>("bsPtMin", 13.0)->setComment("Bs pT must be greater than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);
  }
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1BsMesonSelectionEmulationProducer);
