#include "RecoVertex/HyddraSVProducer/interface/TrackVertexSet.h"

TrackVertexSet::TrackVertexSet(const std::vector<reco::TrackBaseRef> &init, const TransientTrackBuilder* ttBuilder, VertexFitConfig fitConfig) :
  std::set<reco::TrackBaseRef, TrackBaseRefLess>(init.begin(), init.end()),
  ttBuilder_(ttBuilder),
  fitConfig_(fitConfig),
  fitter_(makeFitter(fitConfig_)) {
  fit();
}

TrackVertexSet::TrackVertexSet(std::initializer_list<reco::TrackBaseRef> init, const TransientTrackBuilder* ttBuilder, VertexFitConfig fitConfig) :
  std::set<reco::TrackBaseRef, TrackBaseRefLess>(init),
  ttBuilder_(ttBuilder),
  fitConfig_(fitConfig),
  fitter_(makeFitter(fitConfig_)) {
  fit();
}

TrackVertexSet::TrackVertexSet(const TrackVertexSet& other) :
  std::set<reco::TrackBaseRef, TrackBaseRefLess>(other),
  ttBuilder_(other.ttBuilder_),
  fitConfig_(other.fitConfig_),
  fitter_(makeFitter(fitConfig_))
{fit();}

// Transverse distance between two vertices
double TrackVertexSet::dxySeparation(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparation(this->vertex_, other.vertex_);
}

// Error in transverse distance between two vertices
double TrackVertexSet::dxySeparationError(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparationError(this->vertex_, other.vertex_);
}

// Significance in transverse distance between two vertices
double TrackVertexSet::dxySeparationSignificance(const TrackVertexSet& other) const {
  return VertexHelper::GetDxySeparationSignificance(this->vertex_, other.vertex_);
}

// Distance between two vertices
double TrackVertexSet::distance(const TrackVertexSet& other) const {
  return VertexHelper::GetDistance(this->vertex_, other.vertex_);
}

// Error in the distance between two vertices
double TrackVertexSet::distanceError(const TrackVertexSet& other) const {
  return VertexHelper::GetDistanceError(this->vertex_, other.vertex_);
}

// Significance of the distance between two vertices: distance/distanceError
double TrackVertexSet::distanceSignificance(const TrackVertexSet& other) const {
  return VertexHelper::GetDistanceSignificance(this->vertex_, other.vertex_);
}

// Uses KalmanVertexTrackCompatibilityEstimator to check compatibility between this vertex and one of its tracks
double TrackVertexSet::compatibility(const reco::TrackBaseRef &track) const {

  if(!contains(track)) {
    std::cout << "Warning in TrackVertexSet::compatibility: Given track is not in vertex!" << std::endl;
    return 999.;
  }
  const KalmanVertexTrackCompatibilityEstimator<5> estimator;
  auto result(estimator.estimate(*this, buildTransientTrack(track)));

  if(!result.first)
    throw std::runtime_error("Track compatibility estimation failed");

  return std::sqrt(result.second);
}

double TrackVertexSet::shiftDzAfterTrackRemoval(const reco::TrackBaseRef &track) const {

  TrackVertexSet thisSetMinusTrack = *this;
  thisSetMinusTrack.removeTrack(track);

  return thisSetMinusTrack.isValid()? fabs(this->position().z() - thisSetMinusTrack.position().z()) : -1.;
}

double TrackVertexSet::shift3DAfterTrackRemoval(const reco::TrackBaseRef &track) const {

  TrackVertexSet thisSetMinusTrack = *this;
  thisSetMinusTrack.removeTrack(track);

  return thisSetMinusTrack.isValid()? distance(thisSetMinusTrack) : -1.;
}

std::vector<reco::Track> TrackVertexSet::trackList() const {
  std::vector<reco::TrackBaseRef> returnTracks(this->begin(), this->end());
  std::vector<reco::Track> tracks;
  for(const auto& trackRef : returnTracks)
    tracks.emplace_back(*trackRef);
  return tracks;
}

std::vector<reco::TrackBaseRef> TrackVertexSet::tracks() const {

  std::vector<reco::TrackBaseRef> returnTracks(this->begin(), this->end());
  std::sort(returnTracks.begin(), returnTracks.end(),
    [](const reco::TrackBaseRef& a, const reco::TrackBaseRef& b) {
      return a->pt() < b->pt();
    });
  return returnTracks;
}

// Function to return common tracks between two sets
std::vector<reco::TrackBaseRef> TrackVertexSet::commonTracks(const TrackVertexSet& other) const {
  std::vector<reco::TrackBaseRef> common;
  for (const auto& track : other) {
    if (this->contains(track)) {
      common.emplace_back(track);
    }
  }
  return common;
}

void TrackVertexSet::printTrackInfo() const {

  int count(0);
  std::cout << "\nThis set (normChi2 = " << normChi2() << ") has " << this->size() << " tracks: " << std::endl;
  for (const auto& track : *this) {
    std::cout << "\ttrack " << count << std::endl;
    std::cout << "\t\tpt: " << track->pt() << std::endl;
    std::cout << "\t\teta: " << track->eta() << std::endl;
    std::cout << "\t\tphi: " << track->phi() << std::endl;
    std::cout << "\t\tcompatibility: " << this->compatibility(track) << std::endl;
    count++;
  }
}

// Merge two TrackVertexSets together
TrackVertexSet TrackVertexSet::merge(const TrackVertexSet &other) const {
  return *this + other;
}

// Add a track to the set
void TrackVertexSet::addTrack(const reco::TrackBaseRef& track) {
  this->insert(track);
  fit();
}

// Check if a track exists in the set
bool TrackVertexSet::contains(const reco::TrackBaseRef& track) const {
  return this->find(track) != this->end();
}

// Remove a track from the set
void TrackVertexSet::removeTrack(const reco::TrackBaseRef& track) {
  this->erase(track);
  if(this->size() > 1)
    fit();
  else {
    this->clear();
    vertex_ = TransientVertex();
  }
}

// Clear all tracks from the set
void TrackVertexSet::clearTracks() {
  this->clear();
}

TrackVertexSet& TrackVertexSet::operator=(const TrackVertexSet& other) {
  if (this != &other) {
    std::set<reco::TrackBaseRef, TrackBaseRefLess>::operator=(other);
    ttBuilder_ = other.ttBuilder_;
    fitConfig_ = other.fitConfig_;
    fitter_ = makeFitter(fitConfig_);
    fit();
  }
  return *this;
}

// Overloaded operator to combine two sets into one
TrackVertexSet TrackVertexSet::operator+(const TrackVertexSet& other) const {
  TrackVertexSet result = *this; // Start with the current set
  for(const auto& track : other) {
    result.addTrack(track);
  }

  return result;
}

// Overloaded operator to count overlapping tracks between two sets
size_t TrackVertexSet::operator&(const TrackVertexSet& other) const {
  size_t count = 0;
  for (const auto& track : other) {
    if (this->contains(track)) {
      ++count;
    }
  }
  return count;
}

// Overloaded comparison operator for maintaining uniqueness in a std::set<TrackVertexSet>
bool TrackVertexSet::operator<(const TrackVertexSet& other) const {
    // Compare sizes first
    if (this->size() != other.size())
        return this->size() < other.size();

    // Use the same ref ordering as the underlying set.
    return std::lexicographical_compare(this->begin(), this->end(),
                                        other.begin(), other.end(),
                                        TrackBaseRefLess());
}

// Overloaded equality operator to compare elements
bool TrackVertexSet::operator==(const TrackVertexSet& other) const {
  if(this->size() != other.size()) {
    return false;
  }
  return std::equal(this->begin(), this->end(), other.begin());
}

bool TrackVertexSet::operator!=(const TrackVertexSet& other) const {
  return(!(*this == other));
}

// Overloaded subset operator to check if all tracks in *this are in other
bool TrackVertexSet::operator|=(const TrackVertexSet& other) const {
  if (this->size() > other.size()) {
    return false;
  }
  for (const auto& track : *this) {
    if (!other.contains(track)) {
      return false;
    }
  }
  return true;
}

// private methods
std::unique_ptr<KalmanVertexFitter> TrackVertexSet::makeFitter(const VertexFitConfig& fitConfig) {
  return std::make_unique<KalmanVertexFitter>(fitConfig.useSmoothing, fitConfig.useMuonSystemBounds);
}

void TrackVertexSet::fit() {
  vertex_ = this->size() < 2? TransientVertex() : fitter_->vertex(convertTracks());
}

std::vector<reco::TransientTrack> TrackVertexSet::convertTracks() const {

  std::vector<reco::TransientTrack> ttracks;
  for(const auto& track : *this) {
    ttracks.emplace_back(buildTransientTrack(track));
  }
  return ttracks;
}

reco::TransientTrack TrackVertexSet::buildTransientTrack(const reco::TrackBaseRef& track) const {
  if (dynamic_cast<const reco::GsfTrack*>(&*track))
    return ttBuilder_->build(track.castTo<reco::GsfTrackRef>());
  return ttBuilder_->build(track.castTo<reco::TrackRef>());
}

double TrackVertexSet::calculateChiSquaredPValue(double chiSquaredValue, int degreesOfFreedom) const {

  // Define the chi-squared distribution with the specified degrees of freedom
  boost::math::chi_squared chiSqDist(degreesOfFreedom);

  // Compute the p-value as 1 - CDF(chiSquaredValue)
  double pValue = 1.0 - boost::math::cdf(chiSqDist, chiSquaredValue);

  return pValue;
}

// Direct casting from TrackVertexSet to reco::Vertex
TrackVertexSet::operator reco::Vertex() const {

  if(!vertex_.isValid())
    return reco::Vertex();

  // If smoothing is requested, refit and store
  // vertex-constrained refitted track parameters alongside the original TrackBaseRefs.
  if(fitConfig_.useSmoothing) {
    std::vector<reco::TrackBaseRef> orderedRefs(this->begin(), this->end());
    std::vector<reco::TransientTrack> ttracks;
    ttracks.reserve(orderedRefs.size());
    for(const auto& ref : orderedRefs)
      ttracks.emplace_back(buildTransientTrack(ref));

    KalmanVertexFitter smoother(fitConfig_.useSmoothing, fitConfig_.useMuonSystemBounds);
    TransientVertex smoothed = smoother.vertex(ttracks);

    if(smoothed.isValid() && smoothed.hasRefittedTracks()) {
      reco::Vertex recoVertex(VertexHelper::ConvertFitVertex(smoothed));

      int ifail(0);
      GlobalError err(recoVertex.covariance());
      err.matrix().Inverse(ifail);
      if(!recoVertex.isValid() || ifail != 0)
        return reco::Vertex();

      for(size_t i = 0; i < ttracks.size(); ++i) {
        reco::TransientTrack refitted = smoothed.refittedTrack(ttracks[i]);
        recoVertex.add(orderedRefs[i], refitted.track(), 1.0);
      }
      return recoVertex;
    }
    // Fall through to the unsmoothed path if the smoothed fit fails.
  }

  reco::Vertex recoVertex(VertexHelper::ConvertFitVertex(vertex_));

  // Check for rare edge case where reco::Vertex can fail matrix inversion
  // leading to an exception when checking track compatibility.
  int ifail(0);
  if(this->isValid()) {
    GlobalError err(recoVertex.covariance());
    err.matrix().Inverse(ifail);
  }

  if(!recoVertex.isValid() || ifail != 0)
    return reco::Vertex();

  for(const auto &trackRef : *this) {
    recoVertex.add(trackRef, recoVertex.trackWeight(trackRef));
  }

  return recoVertex;
}
