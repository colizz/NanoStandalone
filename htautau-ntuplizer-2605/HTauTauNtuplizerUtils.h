#ifndef HTAUTAU_NTUPLIZER_UTILS_H
#define HTAUTAU_NTUPLIZER_UTILS_H

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

static constexpr int kMaxObj = 256;
static constexpr float kMissing = -999.f;

// Read a TTreeReaderValue inside const helper functions.
template <typename T>
inline T readerValue(const TTreeReaderValue<T>& v) {
    return *const_cast<TTreeReaderValue<T>&>(v);
}

// Minimal gen-particle record used after copying NanoAOD arrays into a vector.
struct GenParticle {
    float pt = 0, eta = 0, phi = 0, mass = 0;
    int pdgId = 0;
    int status = 0;
    unsigned short statusFlags = 0;
    short mother = -1;
    std::vector<int> dauIdx;

    // NanoAOD statusFlags bit 13 marks the last copy in the generator history.
    bool isLastCopy() const { return statusFlags & (1 << 13); }

    // Build a ROOT four-vector for matching and kinematic calculations.
    TLorentzVector p4() const {
        TLorentzVector v;
        v.SetPtEtaPhiM(pt, eta, phi, mass);
        return v;
    }
};

// Common reconstructed object container for loose leptons and boosted taus.
struct RecoObj {
    int idx = -1;
    float pt = 0, eta = 0, phi = 0, mass = 0;
    int charge = 0;
    // Return the object four-vector using the mass stored in NanoAOD.
    TLorentzVector p4() const {
        TLorentzVector v;
        v.SetPtEtaPhiM(pt, eta, phi, mass);
        return v;
    }
};

// Event-level truth summary for H->tautau and its target boosted channels.
struct GenHTauTau {
    bool hasLHEHiggs = false;
    bool hasHiggs = false;
    bool isSignal = false;
    bool isMuTauH = false;
    bool isETauH = false;
    bool isTauHTauH = false;
    int decayCode = 0; // 0: other/no H, 1: e tauh, 2: mu tauh, 3: tauh tauh
    int tau1Decay = 0; // 1:e, 2:mu, 3:had
    int tau2Decay = 0;
    int tau1VisDecayMode = -1;
    int tau2VisDecayMode = -1;
    GenParticle higgs;
    GenParticle tau1;
    GenParticle tau2;
    GenParticle visTau1;
    GenParticle visTau2;
    GenParticle genEle;
    GenParticle genMu;
    bool hasGenEle = false;
    bool hasGenMu = false;
    bool hasVisTau1 = false;
    bool hasVisTau2 = false;
};

struct GenChannelDR {
    float etauh = kMissing;
    float mutauh = kMissing;
    float tauhtauh = kMissing;
};

// Scalar deltaR implementation from eta/phi coordinates.
inline float deltaR(float eta1, float phi1, float eta2, float phi2) {
    const float deta = eta1 - eta2;
    const float dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    return std::sqrt(deta * deta + dphi * dphi);
}

// Four-vector overload used by reco-gen and reco-reco matching.
inline float deltaR(const TLorentzVector& a, const TLorentzVector& b) {
    return deltaR(a.Eta(), a.Phi(), b.Eta(), b.Phi());
}

// Compute generator-level dR values for the three boosted H->tautau channels.
inline GenChannelDR genChannelDR(const GenHTauTau& gh) {
    GenChannelDR out;
    const GenParticle* visTau = nullptr;
    if (gh.hasVisTau1) visTau = &gh.visTau1;
    if (gh.hasVisTau2) visTau = &gh.visTau2;

    if (gh.isETauH && gh.hasGenEle && visTau) {
        out.etauh = deltaR(gh.genEle.p4(), visTau->p4());
    }
    if (gh.isMuTauH && gh.hasGenMu && visTau) {
        out.mutauh = deltaR(gh.genMu.p4(), visTau->p4());
    }
    if (gh.isTauHTauH && gh.hasVisTau1 && gh.hasVisTau2) {
        out.tauhtauh = deltaR(gh.visTau1.p4(), gh.visTau2.p4());
    }
    return out;
}

// Standard lepton-MET transverse mass used by the e/mu tau_h selections.
inline float transverseMass(float lepPt, float lepPhi, float metPt, float metPhi) {
    const float mt2 = 2.f * lepPt * metPt * (1.f - std::cos(TVector2::Phi_mpi_pi(lepPhi - metPhi)));
    return std::sqrt(std::max(0.f, mt2));
}

// Recursively search a gen-particle decay tree for a descendant with a given abs(pdgId).
inline bool hasDescendant(const std::vector<GenParticle>& gps, int idx, int absPdgId, int* foundIdx = nullptr) {
    if (idx < 0 || idx >= (int)gps.size()) return false;
    for (int dau : gps[idx].dauIdx) {
        if (std::abs(gps[dau].pdgId) == absPdgId) {
            if (foundIdx) *foundIdx = dau;
            return true;
        }
        if (std::abs(gps[dau].pdgId) == std::abs(gps[idx].pdgId)) continue;
        if (hasDescendant(gps, dau, absPdgId, foundIdx)) return true;
    }
    return false;
}

// Classify a generator tau as e, mu, or hadronic using its descendant content.
inline int tauDecayType(const std::vector<GenParticle>& gps, int tauIdx, int* lepIdx = nullptr) {
    int idx = -1;
    if (hasDescendant(gps, tauIdx, 11, &idx)) {
        if (lepIdx) *lepIdx = idx;
        return 1;
    }
    if (hasDescendant(gps, tauIdx, 13, &idx)) {
        if (lepIdx) *lepIdx = idx;
        return 2;
    }
    return 3;
}

// Follow same-PDG daughters until the final copy of a particle is reached.
inline int finalCopyIndex(const std::vector<GenParticle>& gps, int idx) {
    if (idx < 0 || idx >= (int)gps.size()) return idx;
    for (int dau : gps[idx].dauIdx) {
        if (std::abs(gps[dau].pdgId) == std::abs(gps[idx].pdgId)) return finalCopyIndex(gps, dau);
    }
    return idx;
}

// Create the four-vector branch quartet for a common prefix.
inline void setP4Branches(TTree* tree, const std::string& prefix, float& pt, float& eta, float& phi, float& mass) {
    tree->Branch((prefix + "_pt").c_str(), &pt);
    tree->Branch((prefix + "_eta").c_str(), &eta);
    tree->Branch((prefix + "_phi").c_str(), &phi);
    tree->Branch((prefix + "_mass").c_str(), &mass);
}

// Copy a gen-particle four-vector into scalar output variables.
inline void fillP4(const GenParticle& p, float& pt, float& eta, float& phi, float& mass) {
    pt = p.pt; eta = p.eta; phi = p.phi; mass = p.mass;
}

// Copy a reconstructed-object four-vector into scalar output variables.
inline void fillP4(const RecoObj& p, float& pt, float& eta, float& phi, float& mass) {
    pt = p.pt; eta = p.eta; phi = p.phi; mass = p.mass;
}

// Reader bundle for branches common to both ntuplizers.
struct EventBranches {
    TTreeReaderValue<UInt_t> run;
    TTreeReaderValue<UInt_t> luminosityBlock;
    TTreeReaderValue<ULong64_t> event;
    TTreeReaderValue<Float_t> genWeight;
    TTreeReaderValue<Float_t> metPt;
    TTreeReaderValue<Float_t> metPhi;

    // Bind event identifiers and MET branches.
    EventBranches(TTreeReader& r)
        : run(r, "run"), luminosityBlock(r, "luminosityBlock"), event(r, "event"),
          genWeight(r, "genWeight"), metPt(r, "PFMET_pt"), metPhi(r, "PFMET_phi") {}
};

// Reader bundle for truth branches used to build H->tautau labels.
struct GenBranches {
    TTreeReaderValue<Int_t> nGenPart;
    TTreeReaderArray<Float_t> pt, eta, phi, mass;
    TTreeReaderArray<Int_t> pdgId, status;
    TTreeReaderArray<UShort_t> statusFlags;
    TTreeReaderArray<Short_t> mother;

    bool hasLHEPart = false;
    std::unique_ptr<TTreeReaderValue<Int_t>> nLHEPart;
    std::unique_ptr<TTreeReaderArray<Int_t>> lhePdgId;

    TTreeReaderValue<Int_t> nGenVisTau;
    TTreeReaderArray<Float_t> visPt, visEta, visPhi, visMass;
    TTreeReaderArray<Short_t> visMother;
    TTreeReaderArray<UChar_t> visStatus;

    // Bind GenPart and GenVisTau collections; LHEPart is optional in some NanoAOD samples.
    GenBranches(TTreeReader& r)
        : nGenPart(r, "nGenPart"), pt(r, "GenPart_pt"), eta(r, "GenPart_eta"), phi(r, "GenPart_phi"),
          mass(r, "GenPart_mass"), pdgId(r, "GenPart_pdgId"), status(r, "GenPart_status"),
          statusFlags(r, "GenPart_statusFlags"), mother(r, "GenPart_genPartIdxMother"),
          nGenVisTau(r, "nGenVisTau"), visPt(r, "GenVisTau_pt"), visEta(r, "GenVisTau_eta"),
          visPhi(r, "GenVisTau_phi"), visMass(r, "GenVisTau_mass"),
          visMother(r, "GenVisTau_genPartIdxMother"), visStatus(r, "GenVisTau_status") {
        TTree* tree = r.GetTree();
        if (tree && tree->GetBranch("nLHEPart") && tree->GetBranch("LHEPart_pdgId")) {
            hasLHEPart = true;
            nLHEPart = std::make_unique<TTreeReaderValue<Int_t>>(r, "nLHEPart");
            lhePdgId = std::make_unique<TTreeReaderArray<Int_t>>(r, "LHEPart_pdgId");
        }
    }
};

// Build all truth flags and gen four-vectors needed by both reco strategies.
inline GenHTauTau buildGenHTauTau(const GenBranches& br) {
    GenHTauTau out;
    // Keep a separate LHE-level Higgs presence flag when LHEPart is stored in the input.
    if (br.hasLHEPart && br.nLHEPart && br.lhePdgId) {
        for (int i = 0; i < readerValue(*br.nLHEPart); ++i) {
            if (std::abs((*br.lhePdgId)[i]) == 25) out.hasLHEHiggs = true;
        }
    }

    // Copy GenPart arrays into objects and reconstruct mother-daughter links.
    std::vector<GenParticle> gps(readerValue(br.nGenPart));
    for (int i = 0; i < readerValue(br.nGenPart); ++i) {
        gps[i].pt = br.pt[i]; gps[i].eta = br.eta[i]; gps[i].phi = br.phi[i]; gps[i].mass = br.mass[i];
        gps[i].pdgId = br.pdgId[i]; gps[i].status = br.status[i];
        gps[i].statusFlags = br.statusFlags[i]; gps[i].mother = br.mother[i];
    }
    for (int i = 0; i < (int)gps.size(); ++i) {
        if (gps[i].mother >= 0 && gps[i].mother < (int)gps.size()) gps[gps[i].mother].dauIdx.push_back(i);
    }

    // Select the last-copy Higgs as the truth resonance reference.
    int hIdx = -1;
    for (int i = 0; i < (int)gps.size(); ++i) {
        if (std::abs(gps[i].pdgId) == 25 && gps[i].isLastCopy()) {
            hIdx = i;
            break;
        }
    }
    if (hIdx < 0) return out;
    out.hasHiggs = true;
    out.higgs = gps[hIdx];

    // Keep only H->tau tau events for signal-channel labeling.
    std::vector<int> taus;
    for (int dau : gps[hIdx].dauIdx) {
        int fdau = finalCopyIndex(gps, dau);
        if (std::abs(gps[fdau].pdgId) == 15) taus.push_back(fdau);
    }
    if (taus.size() < 2) return out;

    // Determine whether each tau decays to e, mu, or hadrons.
    out.isSignal = true;
    out.tau1 = gps[taus[0]];
    out.tau2 = gps[taus[1]];
    int lepIdx1 = -1, lepIdx2 = -1;
    out.tau1Decay = tauDecayType(gps, taus[0], &lepIdx1);
    out.tau2Decay = tauDecayType(gps, taus[1], &lepIdx2);
    if (out.tau1Decay == 1 && lepIdx1 >= 0) { out.genEle = gps[lepIdx1]; out.hasGenEle = true; }
    if (out.tau2Decay == 1 && lepIdx2 >= 0) { out.genEle = gps[lepIdx2]; out.hasGenEle = true; }
    if (out.tau1Decay == 2 && lepIdx1 >= 0) { out.genMu = gps[lepIdx1]; out.hasGenMu = true; }
    if (out.tau2Decay == 2 && lepIdx2 >= 0) { out.genMu = gps[lepIdx2]; out.hasGenMu = true; }

    // Attach visible hadronic tau information and NanoAOD decay-mode code.
    for (int iv = 0; iv < readerValue(br.nGenVisTau); ++iv) {
        int mother = br.visMother[iv];
        if (mother < 0 || mother >= (int)gps.size()) continue;
        mother = finalCopyIndex(gps, mother);
        GenParticle vis;
        vis.pt = br.visPt[iv]; vis.eta = br.visEta[iv]; vis.phi = br.visPhi[iv]; vis.mass = br.visMass[iv];
        vis.pdgId = gps[mother].pdgId;
        if (mother == taus[0]) {
            out.visTau1 = vis;
            out.tau1VisDecayMode = br.visStatus[iv];
            out.hasVisTau1 = true;
        } else if (mother == taus[1]) {
            out.visTau2 = vis;
            out.tau2VisDecayMode = br.visStatus[iv];
            out.hasVisTau2 = true;
        }
    }

    // Map the tau-pair decay types to the three boosted analysis channels.
    const int a = out.tau1Decay, b = out.tau2Decay;
    out.isETauH = (a == 1 && b == 3) || (a == 3 && b == 1);
    out.isMuTauH = (a == 2 && b == 3) || (a == 3 && b == 2);
    out.isTauHTauH = (a == 3 && b == 3);
    if (out.isETauH) out.decayCode = 1;
    else if (out.isMuTauH) out.decayCode = 2;
    else if (out.isTauHTauH) out.decayCode = 3;
    return out;
}

// Reader bundle for electron and muon branches used in loose lepton selection.
struct LeptonBranches {
    TTreeReaderValue<Int_t> nElectron, nMuon;
    TTreeReaderArray<Float_t> ePt, eEta, ePhi, eMass, eDxy, eDz, eIso;
    TTreeReaderArray<Bool_t> eMvaNoIsoWP90;
    TTreeReaderArray<Int_t> eCharge;
    TTreeReaderArray<Float_t> mPt, mEta, mPhi, mMass, mDxy, mDz, mIso;
    TTreeReaderArray<Bool_t> mLooseId;
    TTreeReaderArray<Int_t> mCharge;

    // Bind only the branches needed by the loose definitions and output p4.
    LeptonBranches(TTreeReader& r)
        : nElectron(r, "nElectron"), nMuon(r, "nMuon"),
          ePt(r, "Electron_pt"), eEta(r, "Electron_eta"), ePhi(r, "Electron_phi"), eMass(r, "Electron_mass"),
          eDxy(r, "Electron_dxy"), eDz(r, "Electron_dz"), eIso(r, "Electron_miniPFRelIso_all"),
          eMvaNoIsoWP90(r, "Electron_mvaNoIso_WP90"), eCharge(r, "Electron_charge"),
          mPt(r, "Muon_pt"), mEta(r, "Muon_eta"), mPhi(r, "Muon_phi"), mMass(r, "Muon_mass"),
          mDxy(r, "Muon_dxy"), mDz(r, "Muon_dz"), mIso(r, "Muon_miniPFRelIso_all"),
          mLooseId(r, "Muon_looseId"), mCharge(r, "Muon_charge") {}
};

// Apply the requested loose electron selection.
inline std::vector<RecoObj> looseElectrons(const LeptonBranches& br) {
    std::vector<RecoObj> out;
    for (int i = 0; i < readerValue(br.nElectron); ++i) {
        if (br.ePt[i] > 10.f && std::abs(br.eEta[i]) < 2.5f && std::abs(br.eDxy[i]) < 0.05f &&
            std::abs(br.eDz[i]) < 0.2f && br.eMvaNoIsoWP90[i] && br.eIso[i] < 0.4f) {
            out.push_back({i, br.ePt[i], br.eEta[i], br.ePhi[i], br.eMass[i], br.eCharge[i]});
        }
    }
    return out;
}

// Apply the requested loose muon selection.
inline std::vector<RecoObj> looseMuons(const LeptonBranches& br) {
    std::vector<RecoObj> out;
    for (int i = 0; i < readerValue(br.nMuon); ++i) {
        if (br.mPt[i] > 10.f && std::abs(br.mEta[i]) < 2.4f && std::abs(br.mDxy[i]) < 0.05f &&
            std::abs(br.mDz[i]) < 0.2f && br.mLooseId[i] && br.mIso[i] < 0.4f) {
            out.push_back({i, br.mPt[i], br.mEta[i], br.mPhi[i], br.mMass[i], br.mCharge[i]});
        }
    }
    return out;
}

// Reader bundle for AK4 jet cleaning, HT, and UParTAK4 b-tag counting.
struct JetBranches {
    TTreeReaderValue<Int_t> nJet;
    TTreeReaderArray<Float_t> pt, eta, phi, mass, rawFactor, btagUParTAK4B;
    TTreeReaderArray<Float_t> chHEF, neHEF, chEmEF, neEmEF;
    TTreeReaderArray<UChar_t> chMultiplicity, neMultiplicity, nConstituents;

    // Bind kinematics, raw factor, UParTAK4 score, and composition-ID inputs.
    JetBranches(TTreeReader& r)
        : nJet(r, "nJet"), pt(r, "Jet_pt"), eta(r, "Jet_eta"), phi(r, "Jet_phi"),
          mass(r, "Jet_mass"), rawFactor(r, "Jet_rawFactor"), btagUParTAK4B(r, "Jet_btagUParTAK4B"),
          chHEF(r, "Jet_chHEF"), neHEF(r, "Jet_neHEF"), chEmEF(r, "Jet_chEmEF"), neEmEF(r, "Jet_neEmEF"),
          chMultiplicity(r, "Jet_chMultiplicity"), neMultiplicity(r, "Jet_neMultiplicity"),
          nConstituents(r, "Jet_nConstituents") {}
};

// Event-level reco quantities shared by both output trees.
struct CommonReco {
    float ht = 0;
    int nLooseEle = 0;
    int nLooseMu = 0;
    int nCleanJet = 0;
    int nBExTight = 0;
    int nBTight = 0;
    int nBMedium = 0;
};

// Build cleaned AK4 jets, HT, and b-tag multiplicities after lepton overlap removal.
inline CommonReco buildCommonReco(const JetBranches& jets, const std::vector<RecoObj>& eles, const std::vector<RecoObj>& mus) {
    CommonReco out;
    out.nLooseEle = eles.size();
    out.nLooseMu = mus.size();
    for (int i = 0; i < readerValue(jets.nJet); ++i) {
        const float rawPt = jets.pt[i] * (1.f - jets.rawFactor[i]);
        if (rawPt <= 25.f || std::abs(jets.eta[i]) >= 2.4f) continue;
        // NanoAOD v15 reference used here does not expose Jet_jetId, so use loose composition cuts.
        const bool passLooseJetId = jets.neHEF[i] < 0.99f && jets.neEmEF[i] < 0.99f &&
                                    jets.nConstituents[i] > 1 && jets.chHEF[i] > 0.f &&
                                    jets.chMultiplicity[i] > 0 && jets.chEmEF[i] < 0.99f;
        if (!passLooseJetId) continue;
        // Remove AK4 jets overlapping loose electrons or muons.
        bool overlaps = false;
        for (const auto& e : eles) overlaps = overlaps || deltaR(jets.eta[i], jets.phi[i], e.eta, e.phi) <= 0.4f;
        for (const auto& m : mus) overlaps = overlaps || deltaR(jets.eta[i], jets.phi[i], m.eta, m.phi) <= 0.4f;
        if (overlaps) continue;
        out.ht += rawPt;
        out.nCleanJet++;
        if (jets.btagUParTAK4B[i] > 0.6298f) out.nBExTight++;
        if (jets.btagUParTAK4B[i] > 0.4648f) out.nBTight++;
        if (jets.btagUParTAK4B[i] > 0.1272f) out.nBMedium++;
    }
    return out;
}

// Add the standard CMS redirector for /store paths while keeping local paths unchanged.
inline std::string inputPathWithRedirector(const char* inputFile) {
    std::string path(inputFile);
    if (path.rfind("/store/", 0) == 0) path = "root://xrootd-cms.infn.it/" + path;
    return path;
}

// Copy run-level generator-weight summary from the NanoAOD Runs tree into output metadata.
inline void writeMetadataTree(TFile* inFile, TFile* outFile) {
    TTree* runs = inFile ? static_cast<TTree*>(inFile->Get("Runs")) : nullptr;
    if (!runs) {
        std::cerr << "Warning: input file has no Runs tree; metadata tree will contain zeros." << std::endl;
    }

    UInt_t run = 0, run_in = 0;
    Long64_t gen_event_count = 0, gen_event_count_in = 0;
    Double_t total_gen_weight = 0., gen_event_sumw_in = 0.;
    Double_t total_gen_weight2 = 0., gen_event_sumw2_in = 0.;

    if (runs) {
        runs->SetBranchAddress("run", &run_in);
        runs->SetBranchAddress("genEventCount", &gen_event_count_in);
        runs->SetBranchAddress("genEventSumw", &gen_event_sumw_in);
        runs->SetBranchAddress("genEventSumw2", &gen_event_sumw2_in);
        for (Long64_t i = 0; i < runs->GetEntries(); ++i) {
            runs->GetEntry(i);
            run = run_in;
            gen_event_count += gen_event_count_in;
            total_gen_weight += gen_event_sumw_in;
            total_gen_weight2 += gen_event_sumw2_in;
        }
    }

    outFile->cd();
    TTree metadata("Metadata", "Run-level metadata copied from NanoAOD Runs tree");
    metadata.Branch("run", &run);
    metadata.Branch("genEventCount", &gen_event_count);
    metadata.Branch("total_gen_weight", &total_gen_weight);
    metadata.Branch("genEventSumw", &total_gen_weight);
    metadata.Branch("genEventSumw2", &total_gen_weight2);
    metadata.Fill();
    metadata.Write();
}

#endif
