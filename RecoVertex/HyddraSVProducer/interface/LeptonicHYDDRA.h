#pragma once

#include <map>
#include <set>
#include <limits>
#include "RecoVertex/HyddraSVProducer/interface/HYDDRABase.h"

// Leptonic displaced vertex reconstruction.
// Targets 2-track vertices from displaced muons and electrons.
// Disambiguation: for each track shared between vertices, keep the vertex
// with the smallest CM cosTheta (most back-to-back decay in the rest frame).
class LeptonicHYDDRA : public HYDDRABase<LeptonicHYDDRA> {

public:

  LeptonicHYDDRA(const edm::ParameterSet& pset) : HYDDRABase(pset) {}

  void seedingImpl(const std::vector<reco::TrackRef>& tracks) {
    this->generateSeeds(tracks);
  }

  void mergingImpl() {
    this->mergeVertices();
  }

  void disambiguationImpl() {
    if (this->empty()) return;

    // Map each track to the vertices that contain it
    std::map<reco::TrackRef, std::vector<const TrackVertexSet*>> trackToVertices;
    for (const auto& vtx : *this)
      for (const auto& track : vtx)
        trackToVertices[track].push_back(&vtx);

    // For each contested track, reject all but the vertex with lowest CM cosTheta
    std::set<const TrackVertexSet*> rejected;

    for (const auto& [track, vtxList] : trackToVertices) {
      std::vector<const TrackVertexSet*> candidates;
      for (const auto* vtx : vtxList)
        if (rejected.find(vtx) == rejected.end()) candidates.push_back(vtx);

      if (candidates.size() < 2) continue;

      const TrackVertexSet* best     = nullptr;
      double                minCosTheta = std::numeric_limits<double>::max();

      for (const auto* vtx : candidates) {
        const double ct = VertexHelper::CalculateCMCosTheta(*vtx, *track);
        if (ct < minCosTheta) { minCosTheta = ct; best = vtx; }
      }

      for (const auto* vtx : candidates)
        if (vtx != best) rejected.insert(vtx);
    }

    if (rejected.empty()) return;

    TrackVertexSetCollection survivors;
    for (const auto& vtx : *this)
      if (rejected.find(&vtx) == rejected.end()) survivors.add(vtx);

    this->clear();
    for (const auto& v : survivors) this->add(v);
  }
};
