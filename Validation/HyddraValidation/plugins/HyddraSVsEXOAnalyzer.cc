// -*- C++ -*-
//
// Package:    RecoVertex/HyddraValidation
// Class:      HyddraSVsEXOAnalyzer
//
// Description: EDAnalyzer for validating HyddraSVsEXOProducer output.
//
//   Reads vertex collections from HyddraSVsEXOProducer and writes a flat TTree.
//   The base collection is selected via the `outputCollection` parameter:
//     "seeds"     — all 2-track seed fits, both pipeline flags stored
//     "inclusive" — tier-0 (seeds -> disambiguation), passDisambiguation trivially 1
//     "isolated"  — tier-1 (seeds -> merge -> disambiguate), both flags trivially 1
//
//   Gen matching (MC only) implements a hybrid approach:
//     Gold:   both tracks match daughters of the same signal gen vertex (ΔR < genDRCut)
//     Bronze: at least one track matches a signal gen daughter
//
// Original Author:  Andres Abreu
//

#include <cmath>
#include <limits>
#include <memory>
#include <vector>

// CMSSW framework
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Data formats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// ROOT
#include "TTree.h"

namespace {

constexpr float INV = -999.f;

float deltaPhi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  while (dphi >  M_PI) dphi -= 2.f * M_PI;
  while (dphi < -M_PI) dphi += 2.f * M_PI;
  return dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {
  float de = eta1 - eta2, dp = deltaPhi(phi1, phi2);
  return std::sqrt(de * de + dp * dp);
}

float computeCosTheta(const reco::Vertex& sv, const reco::Vertex& pv,
                      float px, float py, float pz) {
  float dx = sv.x() - pv.x(), dy = sv.y() - pv.y(), dz = sv.z() - pv.z();
  float disp = std::sqrt(dx*dx + dy*dy + dz*dz);
  float p    = std::sqrt(px*px + py*py + pz*pz);
  if (disp < 1e-9f || p < 1e-9f) return INV;
  return (dx*px + dy*py + dz*pz) / (disp * p);
}

float computeTrackCosThetaCM(float tpx, float tpy, float tpz,
                              float spx, float spy, float spz, float se) {
  if (se < 1e-9f) return INV;
  float bx = spx/se, by = spy/se, bz = spz/se;
  float b2 = bx*bx + by*by + bz*bz;
  if (b2 >= 1.f || b2 < 1e-12f) return INV;
  float gamma = 1.f / std::sqrt(1.f - b2);
  float te    = std::sqrt(tpx*tpx + tpy*tpy + tpz*tpz);
  float bdotp = bx*tpx + by*tpy + bz*tpz;
  float fac   = (gamma - 1.f)*bdotp/b2 - gamma*te;
  float bpx = tpx + fac*bx, bpy = tpy + fac*by, bpz = tpz + fac*bz;
  float bp_mag = std::sqrt(bpx*bpx + bpy*bpy + bpz*bpz);
  float b_mag  = std::sqrt(b2);
  if (bp_mag < 1e-9f) return INV;
  return (bpx*bx + bpy*by + bpz*bz) / (bp_mag * b_mag);
}

float computeDecayAngle(float px1, float py1, float pz1, int q1,
                        float px2, float py2, float pz2, int q2) {
  float spx = px1+px2, spy = py1+py2, spz = pz1+pz2;
  float se  = std::sqrt(px1*px1+py1*py1+pz1*pz1)
            + std::sqrt(px2*px2+py2*py2+pz2*pz2);
  bool use2 = (q2 < q1);
  return computeTrackCosThetaCM(use2?px2:px1, use2?py2:py1, use2?pz2:pz1,
                                 spx, spy, spz, se);
}

float computeDxy(const reco::Vertex& sv, const reco::Vertex& pv) {
  float dx = sv.x()-pv.x(), dy = sv.y()-pv.y();
  return std::sqrt(dx*dx + dy*dy);
}

float computeDxyErr(const reco::Vertex& sv, const reco::Vertex& pv) {
  float dx = sv.x()-pv.x(), dy = sv.y()-pv.y();
  float dxy = std::sqrt(dx*dx + dy*dy);
  if (dxy < 1e-9f) return 0.f;
  float ux = dx/dxy, uy = dy/dxy;
  float sv_err2 = sv.covariance(0,0)*ux*ux + sv.covariance(1,1)*uy*uy + 2.f*sv.covariance(0,1)*ux*uy;
  float pv_err2 = pv.covariance(0,0)*ux*ux + pv.covariance(1,1)*uy*uy + 2.f*pv.covariance(0,1)*ux*uy;
  return std::sqrt(std::max(0.f, sv_err2 + pv_err2));
}

} // namespace

struct HyddraGenVertex {
  float x, y, z, dxy;
  float pt, eta, phi, mass;
  float trk1Pt, trk1Eta, trk1Phi; int trk1PdgId;
  float trk2Pt, trk2Eta, trk2Phi; int trk2PdgId;
  int   motherPdgId;
  std::vector<size_t> dauIdx;
  bool passSelection   = false;
  bool isReconstructed = false;
};

class HyddraSVsEXOAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {

public:
  explicit HyddraSVsEXOAnalyzer(const edm::ParameterSet&);
  ~HyddraSVsEXOAnalyzer() override = default;
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}
  void clearBranches();

  std::vector<HyddraGenVertex> buildGenVertices(const reco::GenParticleCollection&) const;
  std::tuple<int,float,float>  matchTrackToGen(const reco::Track&,
                                                const std::vector<HyddraGenVertex>&,
                                                const reco::GenParticleCollection&) const;
  bool hasSoftMatch(size_t dauIdx, const reco::GenParticleCollection&,
                    const reco::TrackCollection&) const;

  edm::EDGetTokenT<reco::VertexCollection>      seedsToken_;
  edm::EDGetTokenT<reco::VertexCollection>      inclusiveToken_;
  edm::EDGetTokenT<reco::VertexCollection>      isolatedToken_;
  edm::EDGetTokenT<std::vector<int>>            disambFlagsToken_;
  edm::EDGetTokenT<std::vector<int>>            seedIsoFlagsToken_;
  edm::EDGetTokenT<std::vector<int>>            isoFlagsToken_;
  edm::EDGetTokenT<reco::VertexCollection>      pvToken_;
  edm::EDGetTokenT<reco::TrackCollection>       tracksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<pat::METCollection>          metToken_;

  std::string outputCollection_;
  bool  hasGenInfo_;
  int   motherPdgId_;
  float genDRCut_;
  float passSelDRCut_;

  TTree* tree_ = nullptr;

  unsigned int      run_, lumi_;
  unsigned long long event_;
  int   nLeptonTracks_;
  float event_MET_;

  int nHyddraSV_;
  std::vector<float> sv_x_, sv_y_, sv_z_;
  std::vector<float> sv_xErr_, sv_yErr_, sv_zErr_;
  std::vector<float> sv_dxy_, sv_dxyErr_, sv_dxySig_;
  std::vector<float> sv_chi2_, sv_ndof_, sv_normChi2_;
  std::vector<float> sv_pt_, sv_eta_, sv_phi_, sv_mass_, sv_p_;
  std::vector<int>   sv_charge_;
  std::vector<float> sv_cosTheta_, sv_decayAngle_, sv_dR_, sv_beta_;
  std::vector<bool>  sv_passDisambiguation_, sv_passIsolation_;
  std::vector<float> sv_trk1Pt_, sv_trk1Eta_, sv_trk1Phi_;
  std::vector<int>   sv_trk1Charge_;
  std::vector<float> sv_trk1Dxy_, sv_trk1DxyErr_, sv_trk1DxySig_;
  std::vector<float> sv_trk1Dz_, sv_trk1DzErr_, sv_trk1NormChi2_;
  std::vector<float> sv_trk1CosTheta_, sv_trk1CosThetaCM_;
  std::vector<float> sv_trk2Pt_, sv_trk2Eta_, sv_trk2Phi_;
  std::vector<int>   sv_trk2Charge_;
  std::vector<float> sv_trk2Dxy_, sv_trk2DxyErr_, sv_trk2DxySig_;
  std::vector<float> sv_trk2Dz_, sv_trk2DzErr_, sv_trk2NormChi2_;
  std::vector<float> sv_trk2CosTheta_, sv_trk2CosThetaCM_;
  // Gen match (MC only)
  std::vector<float> sv_trk1GenDR_, sv_trk2GenDR_;
  std::vector<float> sv_trk1GenRelPtDiff_, sv_trk2GenRelPtDiff_;
  std::vector<int>   sv_genVtxIdx_, sv_nearestGenVtxIdx_;
  std::vector<float> sv_min3D_;
  std::vector<bool>  sv_isGold_, sv_isBronze_;

  int nHyddraGenSV_;
  std::vector<float> gv_x_, gv_y_, gv_z_, gv_dxy_;
  std::vector<float> gv_pt_, gv_eta_, gv_phi_, gv_mass_;
  std::vector<float> gv_trk1Pt_, gv_trk1Eta_, gv_trk1Phi_; std::vector<int> gv_trk1PdgId_;
  std::vector<float> gv_trk2Pt_, gv_trk2Eta_, gv_trk2Phi_; std::vector<int> gv_trk2PdgId_;
  std::vector<int>   gv_motherPdgId_;
  std::vector<bool>  gv_passSelection_, gv_isReconstructed_;
};

HyddraSVsEXOAnalyzer::HyddraSVsEXOAnalyzer(const edm::ParameterSet& iConfig)
    : seedsToken_       (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("seedVertices")))
    , inclusiveToken_   (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("inclusiveVertices")))
    , isolatedToken_    (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("isolatedVertices")))
    , disambFlagsToken_ (consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("disambiguationFlags")))
    , seedIsoFlagsToken_(consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("seedIsolationFlags")))
    , isoFlagsToken_    (consumes<std::vector<int>>(iConfig.getParameter<edm::InputTag>("isolationFlags")))
    , pvToken_          (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvCollection")))
    , tracksToken_      (consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks")))
    , outputCollection_ (iConfig.getParameter<std::string>("outputCollection"))
    , hasGenInfo_       (iConfig.getParameter<bool>("hasGenInfo"))
    , motherPdgId_      (iConfig.getParameter<int>("motherPdgId"))
    , genDRCut_         (iConfig.getParameter<double>("genDRCut"))
    , passSelDRCut_     (iConfig.getParameter<double>("passSelDRCut"))
{
  usesResource("TFileService");
  if (hasGenInfo_)
    genToken_ = consumes<reco::GenParticleCollection>(
        iConfig.getParameter<edm::InputTag>("genParticles"));
  metToken_ = consumes<pat::METCollection>(
      iConfig.getParameter<edm::InputTag>("MET"));
}

void HyddraSVsEXOAnalyzer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Events", "HyddraSV EXO validation");

  tree_->Branch("run",   &run_);
  tree_->Branch("lumi",  &lumi_);
  tree_->Branch("event", &event_);
  tree_->Branch("nLeptonTracks", &nLeptonTracks_);
  tree_->Branch("Event_MET",     &event_MET_);
  tree_->Branch("nHyddraSV",     &nHyddraSV_);

  tree_->Branch("HyddraSV_x",    &sv_x_);    tree_->Branch("HyddraSV_y",    &sv_y_);    tree_->Branch("HyddraSV_z",    &sv_z_);
  tree_->Branch("HyddraSV_xErr", &sv_xErr_); tree_->Branch("HyddraSV_yErr", &sv_yErr_); tree_->Branch("HyddraSV_zErr", &sv_zErr_);
  tree_->Branch("HyddraSV_dxy",    &sv_dxy_);
  tree_->Branch("HyddraSV_dxyErr", &sv_dxyErr_);
  tree_->Branch("HyddraSV_dxySig", &sv_dxySig_);
  tree_->Branch("HyddraSV_chi2",     &sv_chi2_);
  tree_->Branch("HyddraSV_ndof",     &sv_ndof_);
  tree_->Branch("HyddraSV_normChi2", &sv_normChi2_);
  tree_->Branch("HyddraSV_pt",   &sv_pt_);   tree_->Branch("HyddraSV_eta",  &sv_eta_);
  tree_->Branch("HyddraSV_phi",  &sv_phi_);  tree_->Branch("HyddraSV_mass", &sv_mass_); tree_->Branch("HyddraSV_p", &sv_p_);
  tree_->Branch("HyddraSV_charge",     &sv_charge_);
  tree_->Branch("HyddraSV_cosTheta",   &sv_cosTheta_);
  tree_->Branch("HyddraSV_decayAngle", &sv_decayAngle_);
  tree_->Branch("HyddraSV_dR",         &sv_dR_);
  tree_->Branch("HyddraSV_beta",        &sv_beta_);
  tree_->Branch("HyddraSV_passDisambiguation", &sv_passDisambiguation_);
  tree_->Branch("HyddraSV_passIsolation",      &sv_passIsolation_);
  tree_->Branch("HyddraSV_trk1Pt",         &sv_trk1Pt_);   tree_->Branch("HyddraSV_trk1Eta",    &sv_trk1Eta_);
  tree_->Branch("HyddraSV_trk1Phi",        &sv_trk1Phi_);  tree_->Branch("HyddraSV_trk1Charge", &sv_trk1Charge_);
  tree_->Branch("HyddraSV_trk1Dxy",        &sv_trk1Dxy_);  tree_->Branch("HyddraSV_trk1DxyErr", &sv_trk1DxyErr_);
  tree_->Branch("HyddraSV_trk1DxySig",     &sv_trk1DxySig_);
  tree_->Branch("HyddraSV_trk1Dz",         &sv_trk1Dz_);   tree_->Branch("HyddraSV_trk1DzErr",  &sv_trk1DzErr_);
  tree_->Branch("HyddraSV_trk1NormChi2",   &sv_trk1NormChi2_);
  tree_->Branch("HyddraSV_trk1CosTheta",   &sv_trk1CosTheta_);
  tree_->Branch("HyddraSV_trk1CosThetaCM", &sv_trk1CosThetaCM_);
  tree_->Branch("HyddraSV_trk2Pt",         &sv_trk2Pt_);   tree_->Branch("HyddraSV_trk2Eta",    &sv_trk2Eta_);
  tree_->Branch("HyddraSV_trk2Phi",        &sv_trk2Phi_);  tree_->Branch("HyddraSV_trk2Charge", &sv_trk2Charge_);
  tree_->Branch("HyddraSV_trk2Dxy",        &sv_trk2Dxy_);  tree_->Branch("HyddraSV_trk2DxyErr", &sv_trk2DxyErr_);
  tree_->Branch("HyddraSV_trk2DxySig",     &sv_trk2DxySig_);
  tree_->Branch("HyddraSV_trk2Dz",         &sv_trk2Dz_);   tree_->Branch("HyddraSV_trk2DzErr",  &sv_trk2DzErr_);
  tree_->Branch("HyddraSV_trk2NormChi2",   &sv_trk2NormChi2_);
  tree_->Branch("HyddraSV_trk2CosTheta",   &sv_trk2CosTheta_);
  tree_->Branch("HyddraSV_trk2CosThetaCM", &sv_trk2CosThetaCM_);

  if (hasGenInfo_) {
    tree_->Branch("HyddraSV_trk1GenDR",        &sv_trk1GenDR_);
    tree_->Branch("HyddraSV_trk2GenDR",        &sv_trk2GenDR_);
    tree_->Branch("HyddraSV_trk1GenRelPtDiff", &sv_trk1GenRelPtDiff_);
    tree_->Branch("HyddraSV_trk2GenRelPtDiff", &sv_trk2GenRelPtDiff_);
    tree_->Branch("HyddraSV_genVtxIdx",         &sv_genVtxIdx_);
    tree_->Branch("HyddraSV_nearestGenVtxIdx",  &sv_nearestGenVtxIdx_);
    tree_->Branch("HyddraSV_min3D",             &sv_min3D_);
    tree_->Branch("HyddraSV_isGold",            &sv_isGold_);
    tree_->Branch("HyddraSV_isBronze",          &sv_isBronze_);

    tree_->Branch("nHyddraGenSV",               &nHyddraGenSV_);
    tree_->Branch("HyddraGenSV_x",              &gv_x_);    tree_->Branch("HyddraGenSV_y",  &gv_y_);   tree_->Branch("HyddraGenSV_z",    &gv_z_);
    tree_->Branch("HyddraGenSV_dxy",            &gv_dxy_);
    tree_->Branch("HyddraGenSV_pt",             &gv_pt_);   tree_->Branch("HyddraGenSV_eta", &gv_eta_); tree_->Branch("HyddraGenSV_phi",  &gv_phi_);
    tree_->Branch("HyddraGenSV_mass",           &gv_mass_);
    tree_->Branch("HyddraGenSV_trk1Pt",         &gv_trk1Pt_);  tree_->Branch("HyddraGenSV_trk1Eta", &gv_trk1Eta_);
    tree_->Branch("HyddraGenSV_trk1Phi",        &gv_trk1Phi_); tree_->Branch("HyddraGenSV_trk1PdgId", &gv_trk1PdgId_);
    tree_->Branch("HyddraGenSV_trk2Pt",         &gv_trk2Pt_);  tree_->Branch("HyddraGenSV_trk2Eta", &gv_trk2Eta_);
    tree_->Branch("HyddraGenSV_trk2Phi",        &gv_trk2Phi_); tree_->Branch("HyddraGenSV_trk2PdgId", &gv_trk2PdgId_);
    tree_->Branch("HyddraGenSV_motherPdgId",    &gv_motherPdgId_);
    tree_->Branch("HyddraGenSV_passSelection",  &gv_passSelection_);
    tree_->Branch("HyddraGenSV_isReconstructed",&gv_isReconstructed_);
  }
}

void HyddraSVsEXOAnalyzer::clearBranches() {
  event_MET_ = INV; nLeptonTracks_ = 0;
  sv_x_.clear(); sv_y_.clear(); sv_z_.clear();
  sv_xErr_.clear(); sv_yErr_.clear(); sv_zErr_.clear();
  sv_dxy_.clear(); sv_dxyErr_.clear(); sv_dxySig_.clear();
  sv_chi2_.clear(); sv_ndof_.clear(); sv_normChi2_.clear();
  sv_pt_.clear(); sv_eta_.clear(); sv_phi_.clear(); sv_mass_.clear(); sv_p_.clear();
  sv_charge_.clear(); sv_cosTheta_.clear(); sv_decayAngle_.clear(); sv_dR_.clear(); sv_beta_.clear();
  sv_passDisambiguation_.clear(); sv_passIsolation_.clear();
  sv_trk1Pt_.clear(); sv_trk1Eta_.clear(); sv_trk1Phi_.clear(); sv_trk1Charge_.clear();
  sv_trk1Dxy_.clear(); sv_trk1DxyErr_.clear(); sv_trk1DxySig_.clear();
  sv_trk1Dz_.clear(); sv_trk1DzErr_.clear(); sv_trk1NormChi2_.clear();
  sv_trk1CosTheta_.clear(); sv_trk1CosThetaCM_.clear();
  sv_trk2Pt_.clear(); sv_trk2Eta_.clear(); sv_trk2Phi_.clear(); sv_trk2Charge_.clear();
  sv_trk2Dxy_.clear(); sv_trk2DxyErr_.clear(); sv_trk2DxySig_.clear();
  sv_trk2Dz_.clear(); sv_trk2DzErr_.clear(); sv_trk2NormChi2_.clear();
  sv_trk2CosTheta_.clear(); sv_trk2CosThetaCM_.clear();
  sv_trk1GenDR_.clear(); sv_trk2GenDR_.clear();
  sv_trk1GenRelPtDiff_.clear(); sv_trk2GenRelPtDiff_.clear();
  sv_genVtxIdx_.clear(); sv_nearestGenVtxIdx_.clear(); sv_min3D_.clear();
  sv_isGold_.clear(); sv_isBronze_.clear();
  gv_x_.clear(); gv_y_.clear(); gv_z_.clear(); gv_dxy_.clear();
  gv_pt_.clear(); gv_eta_.clear(); gv_phi_.clear(); gv_mass_.clear();
  gv_trk1Pt_.clear(); gv_trk1Eta_.clear(); gv_trk1Phi_.clear(); gv_trk1PdgId_.clear();
  gv_trk2Pt_.clear(); gv_trk2Eta_.clear(); gv_trk2Phi_.clear(); gv_trk2PdgId_.clear();
  gv_motherPdgId_.clear(); gv_passSelection_.clear(); gv_isReconstructed_.clear();
}

std::vector<HyddraGenVertex>
HyddraSVsEXOAnalyzer::buildGenVertices(const reco::GenParticleCollection& genParts) const {
  std::map<size_t, std::vector<size_t>> motherToDaughters;

  for (size_t i = 0; i < genParts.size(); ++i) {
    const auto& p = genParts[i];
    if (p.status() != 1 || p.charge() == 0) continue;
    const reco::GenParticle* cur = &p;
    for (int step = 0; step < 10; ++step) {
      if (cur->numberOfMothers() == 0) break;
      const reco::GenParticle* mom = dynamic_cast<const reco::GenParticle*>(cur->mother(0));
      if (!mom) break;
      if (std::abs(mom->pdgId()) == motherPdgId_) {
        motherToDaughters[mom - &genParts[0]].push_back(i);
        break;
      }
      if (mom->pdgId() == p.pdgId()) { cur = mom; continue; }
      break;
    }
  }

  std::vector<HyddraGenVertex> result;
  for (auto& [momIdx, dauIdxs] : motherToDaughters) {
    if (dauIdxs.size() < 2) continue;
    float px = 0, py = 0, pz = 0, e = 0;
    for (size_t di : dauIdxs) {
      const auto& d = genParts[di];
      px += d.px(); py += d.py(); pz += d.pz();
      e  += std::sqrt(d.px()*d.px() + d.py()*d.py() + d.pz()*d.pz() + d.mass()*d.mass());
    }
    float p_tot = std::sqrt(px*px + py*py + pz*pz);
    HyddraGenVertex gv;
    gv.x   = genParts[dauIdxs[0]].vx();
    gv.y   = genParts[dauIdxs[0]].vy();
    gv.z   = genParts[dauIdxs[0]].vz();
    gv.dxy = std::sqrt(gv.x*gv.x + gv.y*gv.y);
    gv.pt  = std::sqrt(px*px + py*py);
    gv.eta = (p_tot > 1e-9f) ? std::atanh(pz / p_tot) : 0.f;
    gv.phi = std::atan2(py, px);
    gv.mass = std::sqrt(std::max(0.f, e*e - p_tot*p_tot));
    gv.motherPdgId = genParts[momIdx].pdgId();
    gv.dauIdx = dauIdxs;
    std::sort(gv.dauIdx.begin(), gv.dauIdx.end(),
              [&](size_t a, size_t b){ return genParts[a].pt() > genParts[b].pt(); });
    const auto& d1 = genParts[gv.dauIdx[0]]; const auto& d2 = genParts[gv.dauIdx[1]];
    gv.trk1Pt = d1.pt(); gv.trk1Eta = d1.eta(); gv.trk1Phi = d1.phi(); gv.trk1PdgId = d1.pdgId();
    gv.trk2Pt = d2.pt(); gv.trk2Eta = d2.eta(); gv.trk2Phi = d2.phi(); gv.trk2PdgId = d2.pdgId();
    result.push_back(gv);
  }
  return result;
}

std::tuple<int,float,float>
HyddraSVsEXOAnalyzer::matchTrackToGen(const reco::Track& trk,
                                       const std::vector<HyddraGenVertex>& genVtxs,
                                       const reco::GenParticleCollection& genParts) const {
  int   bestGV  = -1;
  float bestDR  = std::numeric_limits<float>::max();
  float bestRPT = INV;
  for (int gi = 0; gi < (int)genVtxs.size(); ++gi) {
    for (size_t di : genVtxs[gi].dauIdx) {
      const auto& d = genParts[di];
      float dr = deltaR(trk.eta(), trk.phi(), d.eta(), d.phi());
      if (dr < bestDR) {
        bestDR  = dr;
        bestGV  = gi;
        bestRPT = (d.pt() > 1e-9f) ? (trk.pt() - d.pt()) / d.pt() : INV;
      }
    }
  }
  return { bestGV, bestDR, bestRPT };
}

bool HyddraSVsEXOAnalyzer::hasSoftMatch(size_t dauIdx,
                                         const reco::GenParticleCollection& genParts,
                                         const reco::TrackCollection& tracks) const {
  const auto& d = genParts[dauIdx];
  for (const auto& trk : tracks)
    if (deltaR(trk.eta(), trk.phi(), d.eta(), d.phi()) < passSelDRCut_) return true;
  return false;
}

void HyddraSVsEXOAnalyzer::analyze(const edm::Event& iEvent,
                                    const edm::EventSetup&) {
  clearBranches();
  run_   = iEvent.id().run();
  lumi_  = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();

  auto seedsHandle        = iEvent.getHandle(seedsToken_);
  auto inclusiveHandle    = iEvent.getHandle(inclusiveToken_);
  auto isolatedHandle     = iEvent.getHandle(isolatedToken_);
  auto disambFlagsHandle  = iEvent.getHandle(disambFlagsToken_);
  auto seedIsoFlagsHandle = iEvent.getHandle(seedIsoFlagsToken_);
  auto isoFlagsHandle     = iEvent.getHandle(isoFlagsToken_);
  auto pvHandle           = iEvent.getHandle(pvToken_);
  auto tracksHandle       = iEvent.getHandle(tracksToken_);

  nLeptonTracks_ = tracksHandle.isValid() ? static_cast<int>(tracksHandle->size()) : 0;

  edm::Handle<pat::METCollection> metHandle;
  iEvent.getByToken(metToken_, metHandle);
  if (metHandle.isValid() && !metHandle->empty())
    event_MET_ = metHandle->at(0).pt();

  const reco::VertexCollection* svColl     = nullptr;
  const std::vector<int>*       disambFlags = nullptr;
  const std::vector<int>*       isoFlags    = nullptr;
  bool constPassDisamb = false, constPassIso = false;

  if (outputCollection_ == "seeds") {
    svColl      = seedsHandle.isValid()        ? seedsHandle.product()        : nullptr;
    disambFlags = disambFlagsHandle.isValid()  ? disambFlagsHandle.product()  : nullptr;
    isoFlags    = seedIsoFlagsHandle.isValid() ? seedIsoFlagsHandle.product() : nullptr;
  } else if (outputCollection_ == "inclusive") {
    svColl           = inclusiveHandle.isValid() ? inclusiveHandle.product() : nullptr;
    isoFlags         = isoFlagsHandle.isValid()  ? isoFlagsHandle.product()  : nullptr;
    constPassDisamb  = true;
  } else {
    svColl           = isolatedHandle.isValid() ? isolatedHandle.product() : nullptr;
    constPassDisamb  = true;
    constPassIso     = true;
  }

  if (!svColl || pvHandle->empty()) { tree_->Fill(); return; }
  const reco::Vertex& pv = pvHandle->front();

  std::vector<HyddraGenVertex> genVtxs;
  const reco::GenParticleCollection* genParts = nullptr;
  edm::Handle<reco::GenParticleCollection> genHandle;
  if (hasGenInfo_) {
    iEvent.getByToken(genToken_, genHandle);
    if (genHandle.isValid()) {
      genParts = genHandle.product();
      genVtxs  = buildGenVertices(*genParts);
      const reco::TrackCollection& tracks = *tracksHandle;
      for (auto& gv : genVtxs) {
        bool all = true;
        for (size_t di : gv.dauIdx)
          if (!hasSoftMatch(di, *genParts, tracks)) { all = false; break; }
        gv.passSelection = all;
      }
    }
  }

  int svIdx = 0;
  for (const auto& sv : *svColl) {
    std::vector<reco::TrackBaseRef> trkRefs;
    for (auto it = sv.tracks_begin(); it != sv.tracks_end(); ++it)
      trkRefs.push_back(*it);
    if (trkRefs.size() != 2) { ++svIdx; continue; }
    std::sort(trkRefs.begin(), trkRefs.end(),
              [](const reco::TrackBaseRef& a, const reco::TrackBaseRef& b){ return a->pt() > b->pt(); });
    const reco::Track& t1 = *trkRefs[0];
    const reco::Track& t2 = *trkRefs[1];

    bool passDisamb = constPassDisamb ? true
        : (disambFlags && svIdx < (int)disambFlags->size() ? (*disambFlags)[svIdx] != 0 : false);
    bool passIso    = constPassIso    ? true
        : (isoFlags    && svIdx < (int)isoFlags->size()    ? (*isoFlags)[svIdx]    != 0 : false);

    float px = t1.px()+t2.px(), py = t1.py()+t2.py(), pz = t1.pz()+t2.pz();
    float e  = t1.p() +t2.p();
    float p  = std::sqrt(px*px + py*py + pz*pz);
    float pt = std::sqrt(px*px + py*py);
    float eta  = (p > 1e-9f) ? std::atanh(pz/p) : 0.f;
    float phi  = std::atan2(py, px);
    float mass = std::sqrt(std::max(0.f, e*e - p*p));
    float beta_e = std::sqrt(p*p + mass*mass);
    float beta   = (beta_e > 1e-9f) ? p/beta_e : INV;

    float dxy    = computeDxy(sv, pv);
    float dxyErr = computeDxyErr(sv, pv);
    float dxySig = (dxyErr > 1e-9f) ? dxy/dxyErr : 0.f;
    float cosTheta   = computeCosTheta(sv, pv, px, py, pz);
    float decayAngle = computeDecayAngle(t1.px(),t1.py(),t1.pz(),t1.charge(),
                                          t2.px(),t2.py(),t2.pz(),t2.charge());
    float dR     = deltaR(t1.eta(), t1.phi(), t2.eta(), t2.phi());
    int   charge = t1.charge() + t2.charge();

    float t1ct   = (p > 1e-9f) ? (t1.px()*px+t1.py()*py+t1.pz()*pz)/(t1.p()*p) : INV;
    float t2ct   = (p > 1e-9f) ? (t2.px()*px+t2.py()*py+t2.pz()*pz)/(t2.p()*p) : INV;
    float t1ctCM = computeTrackCosThetaCM(t1.px(),t1.py(),t1.pz(), px,py,pz,e);
    float t2ctCM = computeTrackCosThetaCM(t2.px(),t2.py(),t2.pz(), px,py,pz,e);

    float t1dxy = t1.dxy(pv.position()), t1dxyErr = t1.dxyError(pv.position(), pv.covariance());
    float t2dxy = t2.dxy(pv.position()), t2dxyErr = t2.dxyError(pv.position(), pv.covariance());
    float t1dz  = t1.dz(pv.position()),  t1dzErr  = t1.dzError();
    float t2dz  = t2.dz(pv.position()),  t2dzErr  = t2.dzError();

    sv_x_.push_back(sv.x()); sv_y_.push_back(sv.y()); sv_z_.push_back(sv.z());
    sv_xErr_.push_back(sv.xError()); sv_yErr_.push_back(sv.yError()); sv_zErr_.push_back(sv.zError());
    sv_dxy_.push_back(dxy); sv_dxyErr_.push_back(dxyErr); sv_dxySig_.push_back(dxySig);
    sv_chi2_.push_back(sv.chi2()); sv_ndof_.push_back(sv.ndof()); sv_normChi2_.push_back(sv.normalizedChi2());
    sv_pt_.push_back(pt); sv_eta_.push_back(eta); sv_phi_.push_back(phi); sv_mass_.push_back(mass); sv_p_.push_back(p);
    sv_charge_.push_back(charge); sv_cosTheta_.push_back(cosTheta); sv_decayAngle_.push_back(decayAngle);
    sv_dR_.push_back(dR); sv_beta_.push_back(beta);
    sv_passDisambiguation_.push_back(passDisamb); sv_passIsolation_.push_back(passIso);
    sv_trk1Pt_.push_back(t1.pt()); sv_trk1Eta_.push_back(t1.eta()); sv_trk1Phi_.push_back(t1.phi()); sv_trk1Charge_.push_back(t1.charge());
    sv_trk1Dxy_.push_back(t1dxy); sv_trk1DxyErr_.push_back(t1dxyErr); sv_trk1DxySig_.push_back((t1dxyErr>1e-9f)?t1dxy/t1dxyErr:0.f);
    sv_trk1Dz_.push_back(t1dz); sv_trk1DzErr_.push_back(t1dzErr); sv_trk1NormChi2_.push_back(t1.normalizedChi2());
    sv_trk1CosTheta_.push_back(t1ct); sv_trk1CosThetaCM_.push_back(t1ctCM);
    sv_trk2Pt_.push_back(t2.pt()); sv_trk2Eta_.push_back(t2.eta()); sv_trk2Phi_.push_back(t2.phi()); sv_trk2Charge_.push_back(t2.charge());
    sv_trk2Dxy_.push_back(t2dxy); sv_trk2DxyErr_.push_back(t2dxyErr); sv_trk2DxySig_.push_back((t2dxyErr>1e-9f)?t2dxy/t2dxyErr:0.f);
    sv_trk2Dz_.push_back(t2dz); sv_trk2DzErr_.push_back(t2dzErr); sv_trk2NormChi2_.push_back(t2.normalizedChi2());
    sv_trk2CosTheta_.push_back(t2ct); sv_trk2CosThetaCM_.push_back(t2ctCM);

    ++svIdx;

    if (hasGenInfo_ && genParts) {
      auto [gv1, dr1, rpt1] = matchTrackToGen(t1, genVtxs, *genParts);
      auto [gv2, dr2, rpt2] = matchTrackToGen(t2, genVtxs, *genParts);
      sv_trk1GenDR_.push_back(dr1 < 900.f ? dr1 : INV); sv_trk2GenDR_.push_back(dr2 < 900.f ? dr2 : INV);
      sv_trk1GenRelPtDiff_.push_back(rpt1); sv_trk2GenRelPtDiff_.push_back(rpt2);
      bool t1m = (gv1 >= 0 && dr1 < genDRCut_), t2m = (gv2 >= 0 && dr2 < genDRCut_);
      bool isGold = t1m && t2m && (gv1 == gv2);
      if (isGold) genVtxs[gv1].isReconstructed = true;
      sv_isGold_.push_back(isGold); sv_isBronze_.push_back(t1m || t2m);
      sv_genVtxIdx_.push_back(isGold ? gv1 : -1);
      int nearIdx = -1; float minD = std::numeric_limits<float>::max();
      for (int gi = 0; gi < (int)genVtxs.size(); ++gi) {
        float dx=sv.x()-genVtxs[gi].x, dy=sv.y()-genVtxs[gi].y, dz=sv.z()-genVtxs[gi].z;
        float d = std::sqrt(dx*dx+dy*dy+dz*dz);
        if (d < minD) { minD = d; nearIdx = gi; }
      }
      sv_nearestGenVtxIdx_.push_back(nearIdx);
      sv_min3D_.push_back(nearIdx >= 0 ? minD : INV);
    }
  }

  nHyddraSV_ = static_cast<int>(sv_pt_.size());

  if (hasGenInfo_) {
    for (const auto& gv : genVtxs) {
      gv_x_.push_back(gv.x); gv_y_.push_back(gv.y); gv_z_.push_back(gv.z); gv_dxy_.push_back(gv.dxy);
      gv_pt_.push_back(gv.pt); gv_eta_.push_back(gv.eta); gv_phi_.push_back(gv.phi); gv_mass_.push_back(gv.mass);
      gv_trk1Pt_.push_back(gv.trk1Pt); gv_trk1Eta_.push_back(gv.trk1Eta); gv_trk1Phi_.push_back(gv.trk1Phi); gv_trk1PdgId_.push_back(gv.trk1PdgId);
      gv_trk2Pt_.push_back(gv.trk2Pt); gv_trk2Eta_.push_back(gv.trk2Eta); gv_trk2Phi_.push_back(gv.trk2Phi); gv_trk2PdgId_.push_back(gv.trk2PdgId);
      gv_motherPdgId_.push_back(gv.motherPdgId); gv_passSelection_.push_back(gv.passSelection); gv_isReconstructed_.push_back(gv.isReconstructed);
    }
    nHyddraGenSV_ = static_cast<int>(gv_pt_.size());
  }

  tree_->Fill();
}

void HyddraSVsEXOAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("outputCollection", "inclusive");
  desc.add<edm::InputTag>("seedVertices",        edm::InputTag("hyddraSVsEXOProducer", "seedVertices"));
  desc.add<edm::InputTag>("inclusiveVertices",   edm::InputTag("hyddraSVsEXOProducer", "inclusiveVertices"));
  desc.add<edm::InputTag>("isolatedVertices",    edm::InputTag("hyddraSVsEXOProducer", "isolatedVertices"));
  desc.add<edm::InputTag>("disambiguationFlags", edm::InputTag("hyddraSVsEXOProducer", "disambiguationFlags"));
  desc.add<edm::InputTag>("seedIsolationFlags",  edm::InputTag("hyddraSVsEXOProducer", "seedIsolationFlags"));
  desc.add<edm::InputTag>("isolationFlags",      edm::InputTag("hyddraSVsEXOProducer", "isolationFlags"));
  desc.add<edm::InputTag>("pvCollection",        edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("tracks",              edm::InputTag("hyddraSVsEXOProducer", "leptonTracks"));
  desc.add<edm::InputTag>("genParticles",        edm::InputTag("prunedGenParticles"));
  desc.add<edm::InputTag>("MET",                 edm::InputTag("slimmedMETs"));
  desc.add<bool>  ("hasGenInfo",   true);
  desc.add<int>   ("motherPdgId",  9000006);
  desc.add<double>("genDRCut",     0.05);
  desc.add<double>("passSelDRCut", 0.02);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(HyddraSVsEXOAnalyzer);
