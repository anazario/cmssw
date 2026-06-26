// -*- C++ -*-
//
// Validation-only producer that runs HYDDRA directly on PAT electron
// GsfTrackRefs, preserving the GSF edm::Ref identity through TrackBaseRef.

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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoVertex/HyddraSVProducer/interface/LeptonicHYDDRA.h"

class HyddraGsfElectronSVProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraGsfElectronSVProducer(const edm::ParameterSet&);
  ~HyddraGsfElectronSVProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttBuilderToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  LeptonicHYDDRA leptonic_;
};

HyddraGsfElectronSVProducer::HyddraGsfElectronSVProducer(const edm::ParameterSet& iConfig) :
  electronsToken_    (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  pvToken_           (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"))),
  ttBuilderToken_    (esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
  leptonic_          (iConfig.getParameter<edm::ParameterSet>("leptonic"))
{
  produces<reco::TrackCollection> ("leptonTracks");
  produces<reco::VertexCollection>("seedVertices");
  produces<std::vector<int>>      ("disambiguationFlags");
  produces<std::vector<int>>      ("seedIsolationFlags");
  produces<reco::VertexCollection>("inclusiveVertices");
  produces<reco::VertexCollection>("isolatedVertices");
  produces<std::vector<int>>      ("isolationFlags");
}

void HyddraGsfElectronSVProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<pat::ElectronCollection> electronsHandle;
  edm::Handle<reco::VertexCollection> pvHandle;

  iEvent.getByToken(electronsToken_, electronsHandle);
  iEvent.getByToken(pvToken_, pvHandle);

  auto putEmpty = [&]() {
    iEvent.put(std::make_unique<reco::TrackCollection>(),  "leptonTracks");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "seedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "disambiguationFlags");
    iEvent.put(std::make_unique<std::vector<int>>(),       "seedIsolationFlags");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "inclusiveVertices");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "isolatedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "isolationFlags");
  };

  if (!electronsHandle.isValid() || !pvHandle.isValid() || pvHandle->empty()) {
    putEmpty();
    return;
  }

  auto leptonTracks = std::make_unique<reco::TrackCollection>();
  std::vector<reco::TrackBaseRef> trackRefs;
  trackRefs.reserve(electronsHandle->size());

  for (const auto& electron : *electronsHandle) {
    if (electron.gsfTrack().isNull())
      continue;
    trackRefs.emplace_back(electron.gsfTrack());
    leptonTracks->push_back(reco::Track(*electron.gsfTrack()));
  }

  if (trackRefs.size() < 2) {
    iEvent.put(std::move(leptonTracks), "leptonTracks");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "seedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "disambiguationFlags");
    iEvent.put(std::make_unique<std::vector<int>>(),       "seedIsolationFlags");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "inclusiveVertices");
    iEvent.put(std::make_unique<reco::VertexCollection>(), "isolatedVertices");
    iEvent.put(std::make_unique<std::vector<int>>(),       "isolationFlags");
    return;
  }

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(ttBuilderToken_);
  const MagneticField* magneticField = &iSetup.getData(magneticFieldToken_);
  const reco::Vertex& pv = pvHandle->front();

  leptonic_.run_forked(trackRefs, ttBuilder, pv, magneticField);

  iEvent.put(std::move(leptonTracks), "leptonTracks");
  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.seedVertices()),          "seedVertices");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeDisambiguationFlags()), "disambiguationFlags");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeSeedIsolationFlags()),  "seedIsolationFlags");
  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.vertices()),              "inclusiveVertices");
  iEvent.put(std::make_unique<reco::VertexCollection>(leptonic_.isolatedVertices()),      "isolatedVertices");
  iEvent.put(std::make_unique<std::vector<int>>(leptonic_.computeIsolationFlags()),       "isolationFlags");
}

void HyddraGsfElectronSVProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("electrons",    edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("pvCollection", edm::InputTag("offlineSlimmedPrimaryVertices"));

  edm::ParameterSetDescription leptonicDesc;
  leptonicDesc.add<double>("seedCosThetaCut",     -1.0);
  leptonicDesc.add<bool>  ("applySeedChi2Cut",    false);
  leptonicDesc.add<double>("maxNormChi2",          5.0);
  leptonicDesc.add<bool>  ("applyDcaCut",         false);
  leptonicDesc.add<double>("maxDca",              15.0);
  leptonicDesc.add<bool>  ("useSmoothing",        true);
  leptonicDesc.add<bool>  ("useMuonSystemBounds", true);
  desc.add<edm::ParameterSetDescription>("leptonic", leptonicDesc);

  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HyddraGsfElectronSVProducer);
