// -*- C++ -*-
//
// Package:    BsToPhiPhiTo4KFilter
// Class:      BsToPhiPhiTo4KFilter
// 
/**\class BsToPhiPhiTo4KFilter L1Trigger/TTAnalysis/test/BsToPhiPhiTo4KFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
// system include files
#include <memory>
#include "TFile.h"
#include "TH1F.h"

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Gen-level stuff
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Math/interface/deltaPhi.h"

class BsToPhiPhiTo4KFilter: public edm::one::EDFilter<> {
public:

  explicit BsToPhiPhiTo4KFilter(const edm::ParameterSet&);
  ~BsToPhiPhiTo4KFilter() {}
  
private:
  virtual void beginJob() override;
  virtual bool filter(edm::Event&, edm::EventSetup const&) override;
  virtual void endJob() override;

  bool matchParticle(const reco::Candidate* a, const reco::Candidate* b);
  bool passGenFilter(const reco::GenParticleCollection& genParticles, bool verbose=false);
  std::array<unsigned int, 5> evCounters_;
  TH1F *evtCutFlow, *evtCutFlowPhi;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;
  int nEvents_;
};

using std::cout;
using std::endl;
using std::cerr;
using std::setw;

BsToPhiPhiTo4KFilter::BsToPhiPhiTo4KFilter(const edm::ParameterSet& iConfig):
  genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticleInputTag")))
{
}
void BsToPhiPhiTo4KFilter::beginJob() {
  edm::Service<TFileService> fs;

  if (!fs.isAvailable()) return;
  for (size_t i = 0; i < evCounters_.size(); ++i) evCounters_[i] = 0;
  evtCutFlow   = fs->make<TH1F>("evtCutFlow", "",5,-0.5,4.5);
  evtCutFlowPhi   = fs->make<TH1F>("evtCutFlowPhi", "",3,-0.5,2.5);
  nEvents_ = 0;
}
bool BsToPhiPhiTo4KFilter::filter(edm::Event& iEvent, edm::EventSetup const& iSetup) {
  using namespace edm;

  // Apply Gen Filter
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  bool res = iEvent.getByToken(genParticleToken_, genParticleHandle);
  if (res && genParticleHandle.isValid()) {
    reco::GenParticleCollection genParticles = (*genParticleHandle.product());
    ++evCounters_[0];
    evtCutFlow->Fill(0);
    if (!passGenFilter(genParticles, false)) return false;
  }
  else {
    cerr << "filter: GenParticleCollection for InputTag genParticles not found!" << endl; 
    return false;
  }
  ++nEvents_;
  //std::cout<<"Passed Gen"<<std::endl;
  return true;
}
bool BsToPhiPhiTo4KFilter::passGenFilter(const reco::GenParticleCollection& genParticles, bool verbose) {
  auto gbeg = genParticles.begin(); 
  if (verbose) {
    cout << std::setiosflags(std::ios::fixed);
    cout << std::setprecision(2);
    cout << "indx    status    pdgId  charge     eta      phi      pt     energy             mID                             dID"
	 << endl;
  }
  int nphi = 0;  
  size_t i = 0;
  std::array<unsigned int, 3> iCounters = {{0, 0, 0}};
  for (auto it = genParticles.begin(); it != genParticles.end(); ++it) {
    if (std::abs(it->pdgId()) != 333) continue; // Phi
    
    // Find the index of the first mother
    int idm = -1;
    const reco::Candidate* m = it->mother();
    if (m != nullptr) {
      for (auto mit = genParticles.begin(); mit != genParticles.end(); ++mit) {
        const reco::Candidate* ap = &(*mit);
        if (matchParticle(m, ap)) {
	  idm = std::distance(gbeg, mit);
  	  break;
        }
      }
    }
  
    int motherIndex = idm;
    if (motherIndex < 0 || std::abs(m->pdgId()) != 531) continue; // Bs
    
    evtCutFlowPhi->Fill(0);
    std::vector<int> daughterIndices;
    std::ostringstream dID;
    for (size_t j = 0; j < it->numberOfDaughters(); ++j) {
      
      const reco::Candidate* d = it->daughter(j);
      //      if (std::abs(d->pdgId()) != 321 || d->pt() < 2.0 || std::fabs(d->eta()) > 2.4) continue;
      if (std::abs(d->pdgId()) != 321) continue;
      ++iCounters[0];

      if (d->pt() < 2.0) continue;
      ++iCounters[1];

      if (std::fabs(d->eta()) > 2.4) continue;
      ++iCounters[2];
      
      for (auto dit = genParticles.begin(); dit != genParticles.end(); ++dit) {
	const reco::Candidate* ap = &(*dit);  
	if (matchParticle(d, ap)) {
	  int idd = std::distance(gbeg, dit);
	  daughterIndices.push_back(idd);
	  dID << " " << idd;
	  break;
	}
      }
    }
  
    if (daughterIndices.size() < 2) continue;
    evtCutFlowPhi->Fill(1);
    ++nphi;
    if (verbose) {
      std::string ds = dID.str();
      if (!ds.length()) ds = " -";
      cout << setw(4)  << i++
	   << setw(8)  << it->status()
	   << setw(10) << it->pdgId()
	   << setw(8)  << it->charge()
	   << setw(10) << it->eta()
	   << setw(9)  << it->phi()
	   << setw(9)  << it->pt()
	   << setw(9)  << it->energy()
	   << setw(16) << idm
	   << ds
	   << endl;
    }
  }
  cout << std::resetiosflags(std::ios::fixed);

   
  for (size_t i = 0; i < iCounters.size(); ++i) {
    if (iCounters[i] > 3) {
      ++evCounters_[1+i];
      evtCutFlow->Fill(1+i);
    }
  }
  if(nphi>1) {
    ++evCounters_[4];
    evtCutFlow->Fill(4);
    evtCutFlowPhi->Fill(2);
  }
  return (nphi > 1);
}

bool BsToPhiPhiTo4KFilter::matchParticle(const reco::Candidate* a, const reco::Candidate* b) {
  if ( a->pdgId()  == b->pdgId()  &&
       a->status() == b->status() &&        
       std::fabs(a->pt()  - b->pt())  < 1.0e-02 &&        
       std::fabs(a->eta() - b->eta()) < 1.0e-03 &&       
       std::fabs(a->phi() - b->phi()) < 1.0e-03 ) return true;
  return false;
}
void BsToPhiPhiTo4KFilter::endJob() {
  std::vector<std::string> tags {
      "#events ",
      "#phi(mother:Bs) daughters- kaon",
      "#phi(mother:Bs) daughters- kaon pt >=2",
      "#phi(mother:Bs) daughters- kaon eta <=2.4",
      "Passing filter & nPhi > 1",
  };

  std::vector<std::string> tagsPhi {
    "Atleast 1 valid phi with Bs mother",
    "2 prescribed tracks for each phi",
    "Atleast 2 phi",
  };
  /*
  std::cout<<"---------------------------------------------"<<std::endl;
  std::cout<<"------------Cutflow Filter-Evt---------------"<<std::endl;
  std::cout<<"---------------------------------------------"<<std::endl;
  for (size_t i = 0; i < tags.size(); ++i)
    std::cout << setw(40) << tags[i] << setw(10) << evCounters_[i] << std::endl;
  */
  int nbins = evtCutFlow->GetNbinsX();
  for (int i = 1; i <=nbins; ++i) {
    evtCutFlow->GetXaxis()->SetBinLabel(i, tags[i-1].c_str());
  }
  
  int nbinsPhi = evtCutFlowPhi->GetNbinsX();
  for (int i = 1; i <=nbinsPhi; ++i) {
    evtCutFlowPhi->GetXaxis()->SetBinLabel(i, tagsPhi[i-1].c_str());
  }

  std::cout << "BsToPhiPhiTo4KFilter: Events passing filter: " << nEvents_ << std::endl;
}
// define this as a plug-in
DEFINE_FWK_MODULE(BsToPhiPhiTo4KFilter);
