#pragma once

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoVertex/HyddraSVProducer/interface/TrackVertexSetCollection.h"
#include "RecoVertex/HyddraSVProducer/interface/TrackHelper.h"

// CRTP base class for HYDDRA displaced vertex reconstruction.
// Provides shared seeding and merging logic. Derived classes implement
// disambiguationImpl() for their specific vertex selection strategy.
//
// Entry point: run_forked(), which runs two parallel disambiguation paths:
//   Tier 0 (inclusive): seeds -> disambiguation
//   Tier 1 (isolated):  seeds -> merging -> 2-track only -> disambiguation
template <class Derived>
class HYDDRABase : public TrackVertexSetCollection {

protected:

  double seedCosThetaCut_;
  double maxNormChi2_;
  bool applyDcaCut_;
  double maxDca_;
  VertexFitConfig fitConfig_;

  const TransientTrackBuilder* ttBuilder_   = nullptr;
  const MagneticField*         magneticField_ = nullptr;
  const reco::Vertex*          primaryVertex_ = nullptr;

  TrackVertexSetCollection masterList_;
  TrackVertexSetCollection seeds_;
  TrackVertexSetCollection tier0_;
  TrackVertexSetCollection tier1_;

public:

  HYDDRABase(const edm::ParameterSet& pset) {
    seedCosThetaCut_ = pset.getParameter<double>("seedCosThetaCut");
    maxNormChi2_     = pset.getParameter<double>("maxNormChi2");
    applyDcaCut_     = pset.getParameter<bool>("applyDcaCut");
    maxDca_          = pset.getParameter<double>("maxDca");
    fitConfig_.useSmoothing = pset.getParameter<bool>("useSmoothing");
    fitConfig_.useMuonSystemBounds = pset.getParameter<bool>("useMuonSystemBounds");
  }

  // Accessors for the forked pipeline outputs.
  reco::VertexCollection seedVertices()     const { return seeds_.vertices(); }
  reco::VertexCollection isolatedVertices() const { return tier1_.vertices(); }

  // Flag vectors parallel to seedVertices():
  //   [i] = 1 if seed[i] survived tier-0 disambiguation
  std::vector<int> computeDisambiguationFlags() const {
    std::vector<int> flags;
    for (const auto& v : seeds_)
      if (v.isValid()) flags.push_back(tier0_.contains(v) ? 1 : 0);
    return flags;
  }

  // [i] = 1 if seed[i] survived the full isolated (tier-1) path
  std::vector<int> computeSeedIsolationFlags() const {
    std::vector<int> flags;
    for (const auto& v : seeds_)
      if (v.isValid()) flags.push_back(tier1_.contains(v) ? 1 : 0);
    return flags;
  }

  // Flag vector parallel to vertices() (tier-0):
  //   [i] = 1 if inclusive vertex[i] also appears in tier-1
  std::vector<int> computeIsolationFlags() const {
    std::vector<int> flags;
    for (const auto& v : tier0_)
      if (v.isValid()) flags.push_back(tier1_.contains(v) ? 1 : 0);
    return flags;
  }

  // Runs the forked pipeline on the given tracks.
  void run_forked(const std::vector<reco::TrackRef>& tracks,
                  const TransientTrackBuilder* builder,
                  const reco::Vertex& pv,
                  const MagneticField* magneticField) {
    ttBuilder_      = builder;
    magneticField_  = magneticField;
    primaryVertex_  = &pv;
    this->clear();
    masterList_.clear();

    Derived& self = static_cast<Derived&>(*this);

    // Seeding — shared between both tiers
    self.seedingImpl(tracks);
    seeds_ = static_cast<const TrackVertexSetCollection&>(*this);

    // Tier 0: seeds -> disambiguation (inclusive)
    self.disambiguationImpl();
    tier0_ = static_cast<const TrackVertexSetCollection&>(*this);

    // Tier 1: seeds -> merging -> 2-track only -> disambiguation (isolated)
    static_cast<std::set<TrackVertexSet>&>(*this) = seeds_;
    masterList_ = seeds_;
    self.mergingImpl();
    dropMultiTrack();
    self.disambiguationImpl();
    tier1_ = static_cast<const TrackVertexSetCollection&>(*this);

    // Restore *this to tier-0 as the primary output
    static_cast<std::set<TrackVertexSet>&>(*this) = tier0_;
  }

protected:

  // Form all valid 2-track seed vertices from the input tracks.
  void generateSeeds(const std::vector<reco::TrackRef>& tracks) {
    if (tracks.size() < 2) return;

    auto end   = tracks.end();
    auto endm1 = end - 1;

    for (auto x = tracks.begin(); x != endm1; ++x) {
      for (auto y = x + 1; y != end; ++y) {

        if (TrackHelper::OverlappingTrack(**x, **y, ttBuilder_)) continue;

        if (!passesDcaCut(**x, **y)) continue;

        TrackVertexSet seed({*x, *y}, ttBuilder_, fitConfig_);

        if (!isValidVertex(seed)) continue;

        if (primaryVertex_ && primaryVertex_->isValid())
          if (seed.cosTheta(*primaryVertex_) < seedCosThetaCut_) continue;

        this->add(seed);
        masterList_.add(seed);
      }
    }
  }

  // Iteratively merge vertices that share tracks and are spatially close
  // (distance significance < 4), until no further merges are possible.
  void mergeVertices() {
    if (this->empty()) return;

    bool madeChange = true;
    int  iteration  = 0;

    while (madeChange && iteration < 30) {
      madeChange = false;

      TrackVertexSetCollection ignoreList;
      TrackVertexSetCollection mergedList;

      std::vector<TrackVertexSet> vtxVec(this->begin(), this->end());

      for (size_t i = 0; i < vtxVec.size(); ++i) {
        if (ignoreList.contains(vtxVec[i])) continue;

        for (size_t j = i + 1; j < vtxVec.size(); ++j) {
          if (ignoreList.contains(vtxVec[j]) ||
              ignoreList.contains(vtxVec[i])) continue;
          if ((vtxVec[i] & vtxVec[j]) == 0) continue;

          if (vtxVec[i].distanceSignificance(vtxVec[j]) < 4) {
            TrackVertexSet merged(vtxVec[i] + vtxVec[j]);

            if (isValidVertex(merged) && masterList_.doesNotContain(merged)) {
              ignoreList.add(vtxVec[i]);
              ignoreList.add(vtxVec[j]);
              masterList_.add(merged);
              mergedList.add(merged);
              madeChange = true;
            } else if (masterList_.contains(merged)) {
              ignoreList.add(vtxVec[i]);
              ignoreList.add(vtxVec[j]);
            }
          }
        }
      }

      TrackVertexSetCollection updated;
      for (const auto& v : vtxVec)
        if (!ignoreList.contains(v)) updated.add(v);
      updated += mergedList;

      this->clear();
      for (const auto& v : updated) this->add(v);

      ++iteration;
    }
  }

  bool isValidVertex(const TrackVertexSet& set) const {
    return set.isValid() && set.normChi2() < maxNormChi2_;
  }

  bool passesDcaCut(const reco::Track& track1, const reco::Track& track2) const {
    if (!applyDcaCut_) return true;

    TwoTrackMinimumDistance ttmd;
    FreeTrajectoryState fts1(GlobalPoint(track1.vx(), track1.vy(), track1.vz()),
                             GlobalVector(track1.px(), track1.py(), track1.pz()),
                             track1.charge(),
                             magneticField_);
    FreeTrajectoryState fts2(GlobalPoint(track2.vx(), track2.vy(), track2.vz()),
                             GlobalVector(track2.px(), track2.py(), track2.pz()),
                             track2.charge(),
                             magneticField_);
    const bool status = ttmd.calculate(fts1, fts2);
    return !status || ttmd.distance() <= maxDca_;
  }

  // Remove vertices with more than 2 tracks (enforces dilepton constraint in tier-1).
  void dropMultiTrack() {
    TrackVertexSetCollection filtered;
    for (const auto& v : *this)
      if (v.size() == 2) filtered.add(v);
    static_cast<std::set<TrackVertexSet>&>(*this) = filtered;
  }
};
