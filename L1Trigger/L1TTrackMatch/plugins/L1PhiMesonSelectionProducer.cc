// -*- C++ -*-
//
// Package:    L1Trigger/L1TTrackMatch
// Class:      L1PhiMesonSelectionProducer
//
/**\class L1PhiMesonSelectionProducer L1PhiMesonSelectionProducer.cc L1Trigger/L1TTrackMatch/plugins/L1PhiMesonSelectionProducer.cc

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

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidate.h"
#include "DataFormats/L1TCorrelator/interface/TkPhiCandidateFwd.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
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
#include "DataFormats/Math/interface/LorentzVector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TH1D.h>
#include <TFile.h>

//
// class declaration

//
using namespace std;
using namespace edm;
using namespace l1t;

class L1PhiMesonSelectionProducer : public edm::global::EDProducer<> {
public:
  using L1TTTrackType            = TTTrack<Ref_Phase2TrackerDigi_>;
  using TTTrackCollectionType    = std::vector<L1TTTrackType>;
  using TTTrackRef               = edm::Ref<TTTrackCollectionType>;
  using TTTrackRefCollection     = edm::RefVector<TTTrackCollectionType>;
  using TTTrackCollectionHandle  = edm::Handle<TTTrackRefCollection>;
  using TTTrackRefCollectionUPtr = std::unique_ptr<TTTrackRefCollection>;

  explicit L1PhiMesonSelectionProducer(const edm::ParameterSet&);
  ~L1PhiMesonSelectionProducer() override;
  void beginJob() override;
  void endJob();

  void printTrackInfo(const TTTrackRef& track) const;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  static constexpr double KaonMass = 0.493677; // GeV

private:
  // ----------constants, enums and typedefs ---------
  // Relevant constants for the converted track word

  // ----------member functions ----------------------
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  // ----------selectors -----------------------------
  // Based on recommendations from https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGenericSelectors

  // ----------member data ---------------------------
  const edm::EDGetTokenT<TTTrackRefCollection> l1PosKaonTracksToken_;
  const edm::EDGetTokenT<TTTrackRefCollection> l1NegKaonTracksToken_;
  const std::string outputCollectionName_;
  const edm::ParameterSet cutSet_;
  const double dxymax_, dzmax_, dRmax_, tkpairMmin_, tkpairMmax_;
  int debug_;

  TH1D *evtCutFlow, *evtCutFlowTrk;
  mutable std::array<unsigned int, 8> evCounters_;
  mutable std::array<unsigned int, 5> evTrkCounters_;
  TH1D *nPhi;
  TH1D *dz0TrkPairH_;
  TH1D *dRTrkPairH_;
  TH1D *phimassH_;
  TH1D *phimassBinH_;
  TH1D *phimass1H_;
  TH1D *phimass1BinH_;

};

//
// constructors and destructor
//
L1PhiMesonSelectionProducer::L1PhiMesonSelectionProducer(const edm::ParameterSet& iConfig)
  : l1PosKaonTracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1PosKaonTracksInputTag"))),
    l1NegKaonTracksToken_(consumes<TTTrackRefCollection>(iConfig.getParameter<edm::InputTag>("l1NegKaonTracksInputTag"))),
    outputCollectionName_(iConfig.getParameter<std::string>("outputCollectionName")),
                  cutSet_(iConfig.getParameter<edm::ParameterSet>("cutSet")),
                  dxymax_(cutSet_.getParameter<double>("dxymax")),
                  dzmax_(cutSet_.getParameter<double>("dzmax")),
                  dRmax_(cutSet_.getParameter<double>("dRmax")),
              tkpairMmin_(cutSet_.getParameter<double>("tkpairMmin")),
              tkpairMmax_(cutSet_.getParameter<double>("tkpairMmax")),
                   debug_(iConfig.getParameter<int>("debug")) 
{
  produces<TkPhiCandidateCollection>(outputCollectionName_);
  std::cout << "--> L1PhiMesonSelectionProducer: dRmax_ " << dRmax_
            << " dxymax_ " << dxymax_
            << " dzmax_ " << dzmax_
            << " tkpairMmin_ " << tkpairMmin_
            << " tkpairMmax_ " << tkpairMmax_
	    << std::endl;
}

L1PhiMesonSelectionProducer::~L1PhiMesonSelectionProducer() {}
void L1PhiMesonSelectionProducer::beginJob() {
  edm::Service<TFileService> fs;
  if (!fs.isAvailable()) return; 
  evtCutFlow   = fs->make<TH1D>("evtCutFlow", "",8,-0.5,7.5);

  for (size_t i = 0; i < evCounters_.size(); ++i) evCounters_[i] = 0;
  evtCutFlowTrk   = fs->make<TH1D>("evtCutFlowTrk", "",5,-0.5,4.5);
  for (size_t i = 0; i < evTrkCounters_.size(); ++i) evTrkCounters_[i] = 0;

  nPhi      = fs->make<TH1D>("nPhi", "Number of Phi meson candidates", 20, -0.5, 19.5);
  dz0TrkPairH_ = fs->make<TH1D>("dz0TrkPair", "dz0 between Track pairs", 200, -5.0, 5.0);
  dRTrkPairH_  = fs->make<TH1D>("dRTrkPair", "dR between Track pairs", 200, 0.0, 5.0);
  phimassH_    = fs->make<TH1D>("phimass", "Phi Mass", 1000, 0.98, 1.1); //before mass cut
  phimassBinH_    = fs->make<TH1D>("phimassBin", "Phi Mass", 200, 0.98, 1.1);
  phimass1H_    = fs->make<TH1D>("phimass1", "Phi Mass After duplicate removal", 1000, 0.98, 1.1); //after mass cut + dup rem
  phimass1BinH_    = fs->make<TH1D>("phimass1Bin", "Phi Mass dfrnt binning", 200, 0.98, 1.1);
}
void L1PhiMesonSelectionProducer::printTrackInfo(const TTTrackRef& track) const {
  cout << "\t(" 
       << track->momentum().perp() << ", " 
       << track->momentum().eta() << ", " 
       << track->momentum().phi() << ", "
       << track->getStubRefs().size() << ", " 
       << track->stubPtConsistency() << ", " 
       << track->chi2ZRed() << ", "
       << track->chi2XYRed() << ", " 
       << track->z0() << ", " 
       << track->trkMVA1() << ")" 
       << endl;
}
// ------------ method called to produce the data  ------------
void L1PhiMesonSelectionProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const {
  auto L1PhiMesonOutput = std::make_unique<l1t::TkPhiCandidateCollection>();

  TTTrackCollectionHandle l1PosKaonTracksHandle;
  iEvent.getByToken(l1PosKaonTracksToken_, l1PosKaonTracksHandle);

  TTTrackCollectionHandle l1NegKaonTracksHandle;
  iEvent.getByToken(l1NegKaonTracksToken_, l1NegKaonTracksHandle);

  size_t nPosKaonOutputApproximate = l1PosKaonTracksHandle->size();
  size_t nNegKaonOutputApproximate = l1NegKaonTracksHandle->size();
  size_t nPhiMesonOutputApproximate = nPosKaonOutputApproximate + nNegKaonOutputApproximate;
  
  L1PhiMesonOutput->reserve(nPhiMesonOutputApproximate);  

  ++evCounters_[0];
  evtCutFlow->Fill(0);

  if(nPosKaonOutputApproximate > 1 && nNegKaonOutputApproximate > 1) {
    ++evCounters_[1];
    evtCutFlow->Fill(1);
  }
  
  std::array<unsigned int, 4> iCounters = {{0, 0, 0, 0}};
  for (size_t i = 0; i < nPosKaonOutputApproximate; ++i) {
    const auto& trackPosKaonRef = l1PosKaonTracksHandle->at(i);
    const auto& trackPosKaon = *trackPosKaonRef;    
    const edm::Ptr<L1TTTrackType>& trackPosKaonReftoPtr = edm::refToPtr(trackPosKaonRef);

    const GlobalVector& trackPosP = trackPosKaon.momentum();
    math::PtEtaPhiMLorentzVector posKaonP4(trackPosP.perp(), trackPosP.eta(), trackPosP.phi(), KaonMass);

    for (size_t j = 0; j < nNegKaonOutputApproximate; ++j) {
      const auto& trackNegKaonRef = l1NegKaonTracksHandle->at(j);
      const auto& trackNegKaon = *trackNegKaonRef;
      ++evTrkCounters_[0];
      evtCutFlowTrk->Fill(0);
      const edm::Ptr<L1TTTrackType>& trackNegKaonReftoPtr = edm::refToPtr(trackNegKaonRef);      

      const GlobalVector& trackNegP = trackNegKaon.momentum();
      math::PtEtaPhiMLorentzVector negKaonP4(trackNegP.perp(), trackNegP.eta(), trackNegP.phi(), KaonMass);
      
      math::XYZTLorentzVector phiP4(posKaonP4.Px() + negKaonP4.Px(),
				    posKaonP4.Py() + negKaonP4.Py(),
				    posKaonP4.Pz() + negKaonP4.Pz(),
				    posKaonP4.T()  + negKaonP4.T());
    
      TkPhiCandidate tkPhi(phiP4, trackPosKaonReftoPtr, trackNegKaonReftoPtr);

      double dzTrkPair = tkPhi.dzTrkPair();
      dz0TrkPairH_->Fill(dzTrkPair);
      if (std::fabs(dzTrkPair) > dzmax_) continue;
      ++iCounters[0];
      ++evTrkCounters_[1];
      evtCutFlowTrk->Fill(1);

      double dRTrkPair = tkPhi.dRTrkPair();
      dRTrkPairH_->Fill(dRTrkPair);
      if (dRTrkPair > dRmax_) continue;
      ++iCounters[1];
      ++evTrkCounters_[2];
      evtCutFlowTrk->Fill(2);

      
      double mass = tkPhi.p4().M();
      phimassH_->Fill(mass);
      phimassBinH_->Fill(mass);
      
      if (mass < tkpairMmin_ || mass > tkpairMmax_) continue;
      ++iCounters[2];
      ++evTrkCounters_[3];
      evtCutFlowTrk->Fill(3);
      
      bool dupl = false;
      for (const auto& el: *L1PhiMesonOutput) {
        double ptDiff  = el.p4().Pt()  - tkPhi.p4().Pt();
        double etaDiff = el.p4().Eta() - tkPhi.p4().Eta();
        double phiDiff = el.p4().Phi() - tkPhi.p4().Phi();
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
      if (!dupl) {
	++iCounters[3];
	++evTrkCounters_[4];
	evtCutFlowTrk->Fill(4);
	phimass1H_->Fill(mass);
	phimass1BinH_->Fill(mass);
	L1PhiMesonOutput->push_back(tkPhi);      
      }
    }
  }

  for (size_t i = 0; i < iCounters.size(); ++i) {
    if (iCounters[i] > 0) {
      ++evCounters_[2+i];
      evtCutFlow->Fill(2+i);
    }
  }


  if(L1PhiMesonOutput->size() > 0) {
    ++evCounters_[6];
    evtCutFlow->Fill(6);
  }
  if(L1PhiMesonOutput->size() > 1) {
    ++evCounters_[7];
    evtCutFlow->Fill(7);
  }
  nPhi->Fill(L1PhiMesonOutput->size());
  // Put the outputs into the event
  iEvent.put(std::move(L1PhiMesonOutput), outputCollectionName_);
}

void L1PhiMesonSelectionProducer::endJob() {
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

  /*
  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"----------Cutflow Phi-Evt Simulation----------"<<std::endl;
  std::cout<<"---------------------------------------------"<<std::endl;
  for (size_t i = 0; i < tags.size(); ++i)
    std::cout << setw(40) << tags[i] << setw(10) << evCounters_[i] << std::endl;
  /*
  std::cout<<"----------Cutflow Phi-Trk Simulation----------"<<std::endl;
  for (size_t i = 0; i < tagsTrk.size(); ++i)
  std::cout << setw(40) << tagsTrk[i] << setw(10) << evTrkCounters_[i] << std::endl;
  */

  int nbinsPhi = evtCutFlow->GetNbinsX();
  for (int i = 1; i <=nbinsPhi; ++i) {
    evtCutFlow->GetXaxis()->SetBinLabel(i, tags[i-1].c_str());
  }
  
  int nbinsTrk = evtCutFlowTrk->GetNbinsX();
  for (int i = 1; i <=nbinsTrk; ++i) {
    evtCutFlowTrk->GetXaxis()->SetBinLabel(i, tagsTrk[i-1].c_str());
  }
  
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1PhiMesonSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //L1PhiMesonSelectionProducer
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("l1PosKaonTracksInputTag", edm::InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedPositivecharge"));
  desc.add<edm::InputTag>("l1NegKaonTracksInputTag", edm::InputTag("L1KaonTrackSelectionProducer", "Level1TTKaonTracksSelectedNegativecharge"));
  desc.add<std::string>("outputCollectionName", "Level1TTPhiMesonSelected");  {
    edm::ParameterSetDescription descCutSet;
    descCutSet.add<double>("dxymax", 999.0)->setComment("dxy must be less than this value, [cm]");
    descCutSet.add<double>("dzmax", 0.5)->setComment("dz must be less than this value, [cm]");
    descCutSet.add<double>("dRmax", 0.5)->setComment("dr must be less than this value, []");
    descCutSet.add<double>("tkpairMmin", 1.0)->setComment("tkpair mass must be greater than this value, [GeV]");
    descCutSet.add<double>("tkpairMmax", 1.03)->setComment("tkpair mass must be less than this value, [GeV]");
    desc.add<edm::ParameterSetDescription>("cutSet", descCutSet);
  }
  desc.add<int>("debug", 0)->setComment("Verbosity levels: 0, 1, 2, 3");
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1PhiMesonSelectionProducer);
