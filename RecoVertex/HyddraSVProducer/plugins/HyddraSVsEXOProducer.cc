// -*- C++ -*-
//
// Package:    RecoVertex/HyddraSVProducer
// Class:      HyddraSVsEXOProducer
//
// Description: Produces leptonic HYDDRA displaced vertex collections for
//              EXO NanoAOD. Accepts a pre-built reco::TrackCollection (e.g.
//              from HyddraLeptonTrackProducer) and runs the forked pipeline:
//
//   Tier 0 (inclusive): seeds -> disambiguation
//   Tier 1 (isolated):  seeds -> merging -> drop 3+ track -> disambiguation
//
// Original Author:  Andres Abreu
//

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/HyddraSVProducer/interface/LeptonicHYDDRA.h"

class HyddraSVsEXOProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraSVsEXOProducer(const edm::ParameterSet&);
  ~HyddraSVsEXOProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;

  LeptonicHYDDRA leptonic_;
};

HyddraSVsEXOProducer::HyddraSVsEXOProducer(const edm::ParameterSet& iConfig) :
  tracksToken_   (consumes<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("tracks"))),
  pvToken_       (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  leptonic_      (iConfig.getParameter<edm::ParameterSet>("leptonic"))
{
  produces<reco::VertexCollection> ("seedVertices");
  produces<std::vector<int>>       ("disambiguationFlags");
  produces<std::vector<int>>       ("seedIsolationFlags");
  produces<reco::VertexCollection> ("inclusiveVertices");
  produces<reco::VertexCollection> ("isolatedVertices");
  produces<std::vector<int>>       ("isolationFlags");
}

void HyddraSVsEXOProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::TrackCollection> tracksHandle;
  edm::Handle<reco::VertexCollection> pvHandle;

  iEvent.getByToken(tracksToken_, tracksHandle);
  iEvent.getByToken(pvToken_,     pvHandle);

  auto putEmpty = [&]() {
    iEvent.put(std::make_unique<reco::VertexCollection>(), "seedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "disambiguationFlags");
    iEvent.put(std::make_unique<std::vector<int>>(),       "seedIsolationFlags");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "inclusiveVertices");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "isolatedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "isolationFlags");
  };

  if (pvHandle->empty() || tracksHandle->size() < 2) {
    putEmpty();
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
  desc.add<edm::InputTag>("tracks",       edm::InputTag("hyddraLeptonTracks"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlineSlimmedPrimaryVertices"));

  edm::ParameterSetDescription leptonicDesc;
  leptonicDesc.add<double>("seedCosThetaCut", -1.0);
  leptonicDesc.add<double>("maxNormChi2",      5.0);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HyddraSVsEXOProducer);
