// -*- C++ -*-
//
// Package:    RecoVertex/HyddraSVProducer
// Class:      HyddraLeptonTrackProducer
//
// Description: Converts a PAT lepton collection to a reco::TrackCollection
//              for use as input to HyddraSVsEXOProducer.
//
//   leptonType = "muon"     : global -> inner tracker -> tunePMuonBestTrack
//   leptonType = "electron" : GsfTrack sliced to reco::Track
//
// Original Author:  Andres Abreu
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <memory>
#include <string>

class HyddraLeptonTrackProducer : public edm::stream::EDProducer<> {

public:
  explicit HyddraLeptonTrackProducer(const edm::ParameterSet&);
  ~HyddraLeptonTrackProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  enum class Mode { Muon, Electron };

  Mode mode_;
  edm::EDGetTokenT<pat::MuonCollection>     muonsToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
};

HyddraLeptonTrackProducer::HyddraLeptonTrackProducer(const edm::ParameterSet& iConfig) {
  const std::string typeStr = iConfig.getParameter<std::string>("leptonType");
  const edm::InputTag src   = iConfig.getParameter<edm::InputTag>("src");

  if (typeStr == "muon") {
    mode_       = Mode::Muon;
    muonsToken_ = consumes<pat::MuonCollection>(src);
  } else if (typeStr == "electron") {
    mode_          = Mode::Electron;
    electronsToken_ = consumes<pat::ElectronCollection>(src);
  } else {
    throw cms::Exception("Configuration")
        << "HyddraLeptonTrackProducer: leptonType must be 'muon' or 'electron', got '" << typeStr << "'";
  }

  produces<reco::TrackCollection>();
}

void HyddraLeptonTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup&) {
  auto tracks = std::make_unique<reco::TrackCollection>();

  if (mode_ == Mode::Muon) {
    edm::Handle<pat::MuonCollection> muonsHandle;
    iEvent.getByToken(muonsToken_, muonsHandle);
    for (const auto& mu : *muonsHandle) {
      reco::TrackRef tref;
      if      (mu.isGlobalMuon()  && mu.globalTrack().isNonnull()) tref = mu.globalTrack();
      else if (mu.isTrackerMuon() && mu.innerTrack().isNonnull())  tref = mu.innerTrack();
      else if (mu.tunePMuonBestTrack().isNonnull())                tref = mu.tunePMuonBestTrack();
      if (tref.isNonnull()) tracks->push_back(*tref);
    }
  } else {
    edm::Handle<pat::ElectronCollection> electronsHandle;
    iEvent.getByToken(electronsToken_, electronsHandle);
    for (const auto& ele : *electronsHandle) {
      if (ele.gsfTrack().isNonnull())
        tracks->push_back(reco::Track(*ele.gsfTrack()));
    }
  }

  iEvent.put(std::move(tracks));
}

void HyddraLeptonTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>    ("leptonType", "muon");
  desc.add<edm::InputTag>  ("src",        edm::InputTag("slimmedMuons"));
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HyddraLeptonTrackProducer);
