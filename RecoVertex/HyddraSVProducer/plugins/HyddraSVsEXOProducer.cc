// -*- C++ -*-
//
// Package:    RecoVertex/HyddraSVProducer
// Class:      HyddraSVsEXOProducer
//
// Description: Produces leptonic HYDDRA displaced vertex collections for
//              EXO NanoAOD. Accepts PAT muons and electrons, builds a combined
//              reco::TrackCollection, and runs the forked HYDDRA pipeline:
//
//   Tier 0 (inclusive): seeds -> disambiguation
//   Tier 1 (isolated):  seeds -> merging -> drop 3+ track -> disambiguation
//
// Original Author:  Andres Abreu
//

#include <memory>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

// Data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Tracking tools
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// HYDDRA
#include "RecoVertex/HyddraSVProducer/interface/LeptonicHYDDRA.h"

class HyddraSVsEXOProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraSVsEXOProducer(const edm::ParameterSet&);
  ~HyddraSVsEXOProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<pat::MuonCollection>     muonsToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  edm::EDGetTokenT<reco::VertexCollection>  pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  LeptonicHYDDRA leptonic_;
};

HyddraSVsEXOProducer::HyddraSVsEXOProducer(const edm::ParameterSet& iConfig) :
  muonsToken_    (consumes<pat::MuonCollection>    (iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  pvToken_       (consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  leptonic_      (iConfig.getParameter<edm::ParameterSet>("leptonic"))
{
  produces<reco::TrackCollection>  ("leptonTracks");
  produces<reco::VertexCollection> ("seedVertices");
  produces<std::vector<int>>       ("disambiguationFlags");
  produces<std::vector<int>>       ("seedIsolationFlags");
  produces<reco::VertexCollection> ("inclusiveVertices");
  produces<reco::VertexCollection> ("isolatedVertices");
  produces<std::vector<int>>       ("isolationFlags");
}

void HyddraSVsEXOProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<pat::MuonCollection>     muonsHandle;
  edm::Handle<pat::ElectronCollection> electronsHandle;
  edm::Handle<reco::VertexCollection>  pvHandle;

  iEvent.getByToken(muonsToken_,     muonsHandle);
  iEvent.getByToken(electronsToken_, electronsHandle);
  iEvent.getByToken(pvToken_,        pvHandle);

  // Build combined track collection from PAT muons and electrons.
  // Muon track priority: global -> inner tracker -> tunePMuonBestTrack.
  // Electron track: GsfTrack sliced to reco::Track (base class parameters preserved).
  auto tracks = std::make_unique<reco::TrackCollection>();

  for (const auto& mu : *muonsHandle) {
    reco::TrackRef tref;
    if      (mu.isGlobalMuon()  && mu.globalTrack().isNonnull()) tref = mu.globalTrack();
    else if (mu.isTrackerMuon() && mu.innerTrack().isNonnull())  tref = mu.innerTrack();
    else if (mu.tunePMuonBestTrack().isNonnull())                tref = mu.tunePMuonBestTrack();
    if (tref.isNonnull()) tracks->push_back(*tref);
  }

  for (const auto& ele : *electronsHandle) {
    if (ele.gsfTrack().isNonnull())
      tracks->push_back(reco::Track(*ele.gsfTrack()));
  }

  if (tracks->size() > 500)
    edm::LogWarning("HyddraSVsEXOProducer")
      << "Large lepton track collection (" << tracks->size()
      << " tracks). Reconstruction may be slow.";

  // Put tracks first so we can form TrackRefs via the OrphanHandle.
  auto tracksHandle = iEvent.put(std::move(tracks), "leptonTracks");

  // Empty output if no primary vertex or fewer than 2 tracks.
  if (pvHandle->empty() || tracksHandle->size() < 2) {
    iEvent.put(std::make_unique<reco::VertexCollection>(), "seedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "disambiguationFlags");
    iEvent.put(std::make_unique<std::vector<int>>(),       "seedIsolationFlags");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "inclusiveVertices");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "isolatedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "isolationFlags");
    return;
  }

  std::vector<reco::TrackRef> trackRefs;
  trackRefs.reserve(tracksHandle->size());
  for (size_t i = 0; i < tracksHandle->size(); ++i)
    trackRefs.emplace_back(tracksHandle, i);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);
  const reco::Vertex& pv = pvHandle->front();

  leptonic_.run_forked(trackRefs, ttBuilder, pv);

  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.seedVertices()),          "seedVertices");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeDisambiguationFlags()), "disambiguationFlags");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeSeedIsolationFlags()),  "seedIsolationFlags");
  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.vertices()),              "inclusiveVertices");
  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.isolatedVertices()),      "isolatedVertices");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeIsolationFlags()),       "isolationFlags");
}

void HyddraSVsEXOProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<edm::InputTag>("muons",        edm::InputTag("slimmedMuons"));
  desc.add<edm::InputTag>("electrons",    edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlineSlimmedPrimaryVertices"));

  edm::ParameterSetDescription leptonicDesc;
  leptonicDesc.add<double>("seedCosThetaCut", -1.0);
  leptonicDesc.add<double>("maxNormChi2",      5.0);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HyddraSVsEXOProducer);
