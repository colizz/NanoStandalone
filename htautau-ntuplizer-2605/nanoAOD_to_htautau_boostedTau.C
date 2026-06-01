/*
 * BoostedTau H->tautau ntuplizer.
 *
 * Purpose:
 *   Build a compact event-level ntuple for the boostedTau-based H->tautau
 *   reconstruction strategy. The script evaluates three mutually relevant
 *   final-state hypotheses in each selected event: mu tau_h, e tau_h, and
 *   tau_h tau_h. For each channel it stores the best passing candidate, plus
 *   common event, reco-summary, generator-truth, and metadata information.
 *
 * Input and output:
 *   - Reads NanoAOD Events and Runs trees. Plain /store paths are resolved
 *     through the XRootD redirector in the shared utility header.
 *   - Writes an Events-like tree named HTauTauBoostedTau.
 *   - Writes a Metadata tree from Runs with genEventCount, total gen weight
 *     / genEventSumw, and genEventSumw2 when available.
 *   - Per-event genWeight from Events is stored in the main ntuple.
 *
 * Event and common object selections:
 *   - Events must contain at least one boostedTau with pT > 30 GeV and
 *     |eta| < 2.4; otherwise the event is skipped.
 *   - Loose electrons: pt > 10, |eta| < 2.5, |dxy| < 0.05, |dz| < 0.2,
 *     mvaNoIso_WP90, miniPFRelIso_all < 0.4.
 *   - Loose muons: pt > 10, |eta| < 2.4, |dxy| < 0.05, |dz| < 0.2,
 *     looseId, miniPFRelIso_all < 0.4.
 *   - Cleaned AK4 jets: raw pT > 25 GeV, |eta| < 2.4, loose jet ID, and
 *     dR > 0.4 from all loose electrons and loose muons.
 *   - HT is the scalar pT sum of cleaned AK4 jets. b-tag counts use
 *     Jet_btagUParTAK4B with WPs: medium 0.1272, tight 0.4648,
 *     ex-tight 0.6298.
 *
 * Channel selections:
 *   - mu tau_h: require no loose electrons and at least one loose muon.
 *     Candidate pairs require 0.1 < dR(mu, tau) < 0.8, opposite sign,
 *     MT(mu, MET) < 80 GeV, and vector-sum Higgs pT(mu + tau + MET) > 200 GeV.
 *   - e tau_h: require no loose muons and at least one loose electron.
 *     Candidate pairs require 0.1 < dR(e, tau) < 0.8, opposite sign,
 *     MT(e, MET) < 80 GeV, and vector-sum Higgs pT(e + tau + MET) > 200 GeV.
 *   - tau_h tau_h: require no loose electrons, no loose muons, and at least
 *     two selected boostedTau objects. Candidate pairs require
 *     0.1 < dR(tau1, tau2) < 0.8, opposite sign, and vector-sum
 *     Higgs pT(tau1 + tau2 + MET) > 200 GeV.
 *
 * Best-candidate rule:
 *   object_selection_rule is exposed as the third macro argument.
 *   - "max_pt_higgs" (default): choose the passing pair with largest
 *     reconstructed Higgs pT in each channel.
 *   - "gen_matched": prefer candidates matched to the corresponding gen
 *     e/mu and visible hadronic tau(s); if several candidates are matched,
 *     choose the one with largest reconstructed Higgs pT. If none are matched,
 *     fall back to largest reconstructed Higgs pT.
 *
 * Generator truth and matching:
 *   - Finds the last-copy gen Higgs from GenPart and labels H->tautau signal.
 *   - Tau decays are classified as e tau_h, mu tau_h, tau_h tau_h, or other.
 *   - LHEPart is optional: has_lhe_higgs is filled only when LHE branches exist.
 *   - gen_etauh_dr, gen_mutauh_dr, and gen_tauhtauh_dr are filled only for the
 *     corresponding gen channel using gen e/mu and visible hadronic tau objects;
 *     otherwise they are set to kMissing.
 *   - Channel gen-match flags are evaluated only for signal events in the
 *     matching gen channel, requiring dR < 0.4 for the relevant reco objects.
 *
 * Stored variable groups:
 *   - Event: run, luminosityBlock, event, genWeight, MET_pt, MET_phi.
 *   - Reco summary: HT, loose lepton counts, boostedTau count, cleaned AK4 jet
 *     count, cleaned AK4 b-tag counts at medium/tight/ex-tight WPs.
 *   - Gen: Higgs/tau four-vectors, decay labels, tau decay modes, visible tau
 *     decay modes, channel dR values, and signal/channel flags.
 *   - mu tau_h: reco flag, gen-match flag, Higgs pT, pair dR, MT, muon
 *     four-vector and charge, and all stored NanoAOD boostedTau variables for
 *     the selected tau.
 *   - e tau_h: reco flag, gen-match flag, Higgs pT, pair dR, MT, electron
 *     four-vector and charge, and all stored NanoAOD boostedTau variables for
 *     the selected tau.
 *   - tau_h tau_h: reco flag, gen-match flag, Higgs pT, pair dR, and all
 *     stored NanoAOD boostedTau variables for both selected taus.
 *
 * Notes:
 *   - Non-selected channels are kept in the row with reco/gen-match flags set
 *     to 0 and kinematic/object branches filled with dummy values.
 *   - Dummy float values use kMissing from HTauTauNtuplizerUtils.h.
 */
#include "HTauTauNtuplizerUtils.h"

#include "TFile.h"

// Reader bundle for all boostedTau branches copied to the channel outputs.
struct BoostedTauBranches {
    TTreeReaderValue<Int_t> n;
    TTreeReaderArray<Float_t> pt, eta, phi, mass;
    TTreeReaderArray<Int_t> charge, decayMode;
    TTreeReaderArray<UChar_t> genPartFlav, idAntiEle2018, idAntiMu, idMVAnewDM2017v2, idMVAoldDM2017v2;
    TTreeReaderArray<Short_t> genPartIdx, jetIdx, rawAntiEleCat2018;
    TTreeReaderArray<Float_t> chargedIso, neutralIso, leadTkDeltaEta, leadTkDeltaPhi, leadTkPtOverTauPt;
    TTreeReaderArray<Float_t> photonsOutsideSignalCone, puCorr, rawAntiEle2018;
    TTreeReaderArray<Float_t> rawBoostedDeepTauRunIIv2p0VSe, rawBoostedDeepTauRunIIv2p0VSjet, rawBoostedDeepTauRunIIv2p0VSmu;
    TTreeReaderArray<Float_t> rawIso, rawIsodR03, rawMVAnewDM2017v2, rawMVAoldDM2017v2;

    // Bind boostedTau kinematics, IDs, isolation, gen matching, and raw discriminator inputs.
    BoostedTauBranches(TTreeReader& r)
        : n(r, "nboostedTau"), pt(r, "boostedTau_pt"), eta(r, "boostedTau_eta"), phi(r, "boostedTau_phi"),
          mass(r, "boostedTau_mass"), charge(r, "boostedTau_charge"), decayMode(r, "boostedTau_decayMode"),
          genPartFlav(r, "boostedTau_genPartFlav"), idAntiEle2018(r, "boostedTau_idAntiEle2018"),
          idAntiMu(r, "boostedTau_idAntiMu"), idMVAnewDM2017v2(r, "boostedTau_idMVAnewDM2017v2"),
          idMVAoldDM2017v2(r, "boostedTau_idMVAoldDM2017v2"), genPartIdx(r, "boostedTau_genPartIdx"),
          jetIdx(r, "boostedTau_jetIdx"), rawAntiEleCat2018(r, "boostedTau_rawAntiEleCat2018"),
          chargedIso(r, "boostedTau_chargedIso"), neutralIso(r, "boostedTau_neutralIso"),
          leadTkDeltaEta(r, "boostedTau_leadTkDeltaEta"), leadTkDeltaPhi(r, "boostedTau_leadTkDeltaPhi"),
          leadTkPtOverTauPt(r, "boostedTau_leadTkPtOverTauPt"),
          photonsOutsideSignalCone(r, "boostedTau_photonsOutsideSignalCone"), puCorr(r, "boostedTau_puCorr"),
          rawAntiEle2018(r, "boostedTau_rawAntiEle2018"),
          rawBoostedDeepTauRunIIv2p0VSe(r, "boostedTau_rawBoostedDeepTauRunIIv2p0VSe"),
          rawBoostedDeepTauRunIIv2p0VSjet(r, "boostedTau_rawBoostedDeepTauRunIIv2p0VSjet"),
          rawBoostedDeepTauRunIIv2p0VSmu(r, "boostedTau_rawBoostedDeepTauRunIIv2p0VSmu"),
          rawIso(r, "boostedTau_rawIso"), rawIsodR03(r, "boostedTau_rawIsodR03"),
          rawMVAnewDM2017v2(r, "boostedTau_rawMVAnewDM2017v2"),
          rawMVAoldDM2017v2(r, "boostedTau_rawMVAoldDM2017v2") {}
};

// Output buffer for one selected boostedTau object.
struct TauOut {
    Float_t pt, eta, phi, mass;
    Int_t charge, decayMode, genPartFlav;
    Int_t genPartIdx, jetIdx, idAntiEle2018, idAntiMu, idMVAnewDM2017v2, idMVAoldDM2017v2;
    Int_t rawAntiEleCat2018;
    Float_t chargedIso, neutralIso, leadTkDeltaEta, leadTkDeltaPhi, leadTkPtOverTauPt;
    Float_t photonsOutsideSignalCone, puCorr, rawAntiEle2018;
    Float_t rawBoostedDeepTauRunIIv2p0VSe, rawBoostedDeepTauRunIIv2p0VSjet, rawBoostedDeepTauRunIIv2p0VSmu;
    Float_t rawIso, rawIsodR03, rawMVAnewDM2017v2, rawMVAoldDM2017v2;

    // Reset the boostedTau output block when a channel has no selected candidate.
    void reset() {
        pt = eta = phi = mass = kMissing;
        charge = decayMode = genPartFlav = genPartIdx = jetIdx = idAntiEle2018 = idAntiMu = idMVAnewDM2017v2 = idMVAoldDM2017v2 = rawAntiEleCat2018 = -999;
        chargedIso = neutralIso = leadTkDeltaEta = leadTkDeltaPhi = leadTkPtOverTauPt = kMissing;
        photonsOutsideSignalCone = puCorr = rawAntiEle2018 = kMissing;
        rawBoostedDeepTauRunIIv2p0VSe = rawBoostedDeepTauRunIIv2p0VSjet = rawBoostedDeepTauRunIIv2p0VSmu = kMissing;
        rawIso = rawIsodR03 = rawMVAnewDM2017v2 = rawMVAoldDM2017v2 = kMissing;
    }
};

// Register the full boostedTau output block with a channel-specific prefix.
void branchTau(TTree* t, const std::string& p, TauOut& o) {
    t->Branch((p + "_pt").c_str(), &o.pt);
    t->Branch((p + "_eta").c_str(), &o.eta);
    t->Branch((p + "_phi").c_str(), &o.phi);
    t->Branch((p + "_mass").c_str(), &o.mass);
    t->Branch((p + "_charge").c_str(), &o.charge);
    t->Branch((p + "_decayMode").c_str(), &o.decayMode);
    t->Branch((p + "_genPartFlav").c_str(), &o.genPartFlav);
    t->Branch((p + "_genPartIdx").c_str(), &o.genPartIdx);
    t->Branch((p + "_jetIdx").c_str(), &o.jetIdx);
    t->Branch((p + "_idAntiEle2018").c_str(), &o.idAntiEle2018);
    t->Branch((p + "_idAntiMu").c_str(), &o.idAntiMu);
    t->Branch((p + "_idMVAnewDM2017v2").c_str(), &o.idMVAnewDM2017v2);
    t->Branch((p + "_idMVAoldDM2017v2").c_str(), &o.idMVAoldDM2017v2);
    t->Branch((p + "_chargedIso").c_str(), &o.chargedIso);
    t->Branch((p + "_neutralIso").c_str(), &o.neutralIso);
    t->Branch((p + "_leadTkDeltaEta").c_str(), &o.leadTkDeltaEta);
    t->Branch((p + "_leadTkDeltaPhi").c_str(), &o.leadTkDeltaPhi);
    t->Branch((p + "_leadTkPtOverTauPt").c_str(), &o.leadTkPtOverTauPt);
    t->Branch((p + "_photonsOutsideSignalCone").c_str(), &o.photonsOutsideSignalCone);
    t->Branch((p + "_puCorr").c_str(), &o.puCorr);
    t->Branch((p + "_rawAntiEle2018").c_str(), &o.rawAntiEle2018);
    t->Branch((p + "_rawAntiEleCat2018").c_str(), &o.rawAntiEleCat2018);
    t->Branch((p + "_rawBoostedDeepTauRunIIv2p0VSe").c_str(), &o.rawBoostedDeepTauRunIIv2p0VSe);
    t->Branch((p + "_rawBoostedDeepTauRunIIv2p0VSjet").c_str(), &o.rawBoostedDeepTauRunIIv2p0VSjet);
    t->Branch((p + "_rawBoostedDeepTauRunIIv2p0VSmu").c_str(), &o.rawBoostedDeepTauRunIIv2p0VSmu);
    t->Branch((p + "_rawIso").c_str(), &o.rawIso);
    t->Branch((p + "_rawIsodR03").c_str(), &o.rawIsodR03);
    t->Branch((p + "_rawMVAnewDM2017v2").c_str(), &o.rawMVAnewDM2017v2);
    t->Branch((p + "_rawMVAoldDM2017v2").c_str(), &o.rawMVAoldDM2017v2);
}

// Copy one NanoAOD boostedTau entry into the output buffer.
void fillTau(const BoostedTauBranches& br, int idx, TauOut& o) {
    o.pt = br.pt[idx]; o.eta = br.eta[idx]; o.phi = br.phi[idx]; o.mass = br.mass[idx];
    o.charge = br.charge[idx]; o.decayMode = br.decayMode[idx]; o.genPartFlav = br.genPartFlav[idx];
    o.genPartIdx = br.genPartIdx[idx]; o.jetIdx = br.jetIdx[idx]; o.idAntiEle2018 = br.idAntiEle2018[idx];
    o.idAntiMu = br.idAntiMu[idx]; o.idMVAnewDM2017v2 = br.idMVAnewDM2017v2[idx];
    o.idMVAoldDM2017v2 = br.idMVAoldDM2017v2[idx]; o.chargedIso = br.chargedIso[idx];
    o.neutralIso = br.neutralIso[idx]; o.leadTkDeltaEta = br.leadTkDeltaEta[idx];
    o.leadTkDeltaPhi = br.leadTkDeltaPhi[idx]; o.leadTkPtOverTauPt = br.leadTkPtOverTauPt[idx];
    o.photonsOutsideSignalCone = br.photonsOutsideSignalCone[idx]; o.puCorr = br.puCorr[idx];
    o.rawAntiEle2018 = br.rawAntiEle2018[idx]; o.rawAntiEleCat2018 = br.rawAntiEleCat2018[idx];
    o.rawBoostedDeepTauRunIIv2p0VSe = br.rawBoostedDeepTauRunIIv2p0VSe[idx];
    o.rawBoostedDeepTauRunIIv2p0VSjet = br.rawBoostedDeepTauRunIIv2p0VSjet[idx];
    o.rawBoostedDeepTauRunIIv2p0VSmu = br.rawBoostedDeepTauRunIIv2p0VSmu[idx];
    o.rawIso = br.rawIso[idx]; o.rawIsodR03 = br.rawIsodR03[idx];
    o.rawMVAnewDM2017v2 = br.rawMVAnewDM2017v2[idx]; o.rawMVAoldDM2017v2 = br.rawMVAoldDM2017v2[idx];
}

// Match a reconstructed boosted tau to either visible hadronic gen tau.
bool tauMatchesAnyVisible(const RecoObj& tau, const GenHTauTau& gh) {
    bool ok = false;
    if (gh.hasVisTau1) ok = ok || deltaR(tau.p4(), gh.visTau1.p4()) < 0.4f;
    if (gh.hasVisTau2) ok = ok || deltaR(tau.p4(), gh.visTau2.p4()) < 0.4f;
    return ok;
}

// Compare two passing candidates according to the requested best-object rule.
bool preferCandidate(const std::string& rule, float candHpt, bool candGenMatched, float bestHpt, bool bestGenMatched) {
    if (rule == "gen_matched") {
        if (candGenMatched != bestGenMatched) return candGenMatched;
    }
    return candHpt > bestHpt;
}

// Convert NanoAOD events into an event-level ntuple for boosted tau reconstruction.
void nanoAOD_to_htautau_boostedTau(const char* inputFile, const char* outputFile, const char* object_selection_rule = "max_pt_higgs") {
    const std::string selectionRule(object_selection_rule);
    if (selectionRule != "max_pt_higgs" && selectionRule != "gen_matched") {
        std::cerr << "Unsupported object_selection_rule=\"" << selectionRule
                  << "\". Supported values are max_pt_higgs and gen_matched." << std::endl;
        return;
    }

    // Resolve /store paths through XRootD, then open the NanoAOD input.
    const std::string inputPath = inputPathWithRedirector(inputFile);
    TFile* inFile = TFile::Open(inputPath.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return;
    }

    // Bind input branches for event info, truth, loose leptons, AK4 jets, and boosted taus.
    TTreeReader reader("Events", inFile);
    EventBranches ev(reader);
    GenBranches gen(reader);
    LeptonBranches lep(reader);
    JetBranches jets(reader);
    BoostedTauBranches bt(reader);

    // Create one event-level tree containing three independent boosted-tau channels.
    TFile* outFile = new TFile(outputFile, "RECREATE");
    TTree* out = new TTree("HTauTauBoostedTau", "H->tautau boosted tau ntuple");

    // Output storage: common event variables, truth labels, and per-channel candidates.
    UInt_t out_run = 0, out_lumi = 0;
    ULong64_t out_event = 0;
    Float_t gen_weight = 0, met_pt = 0, met_phi = 0, ht = 0;
    Int_t n_loose_e = 0, n_loose_mu = 0, n_boosted_tau = 0, n_clean_jet = 0;
    Int_t n_b_extight = 0, n_b_tight = 0, n_b_medium = 0;
    Int_t has_lhe_h = 0, has_gen_h = 0, is_sig = 0, is_mutauh = 0, is_etauh = 0, is_tauhtauh = 0, gen_decay = 0;
    Float_t gen_h_pt = kMissing, gen_h_eta = kMissing, gen_h_phi = kMissing, gen_h_mass = kMissing;
    Float_t gen_tau1_pt = kMissing, gen_tau1_eta = kMissing, gen_tau1_phi = kMissing, gen_tau1_mass = kMissing;
    Float_t gen_tau2_pt = kMissing, gen_tau2_eta = kMissing, gen_tau2_phi = kMissing, gen_tau2_mass = kMissing;
    Float_t gen_etauh_dr = kMissing, gen_mutauh_dr = kMissing, gen_tauhtauh_dr = kMissing;
    Int_t gen_tau1_decay = 0, gen_tau2_decay = 0, gen_tau1_vis_dm = -1, gen_tau2_vis_dm = -1;
    Int_t mutauh_reco = 0, mutauh_gen_match = 0, etauh_reco = 0, etauh_gen_match = 0, tauhtauh_reco = 0, tauhtauh_gen_match = 0;
    Float_t mutauh_higgs_pt = kMissing, mutauh_dr = kMissing, mutauh_mt = kMissing;
    Float_t etauh_higgs_pt = kMissing, etauh_dr = kMissing, etauh_mt = kMissing;
    Float_t tauhtauh_higgs_pt = kMissing, tauhtauh_dr = kMissing;
    Float_t mu_pt = kMissing, mu_eta = kMissing, mu_phi = kMissing, mu_mass = kMissing;
    Float_t ele_pt = kMissing, ele_eta = kMissing, ele_phi = kMissing, ele_mass = kMissing;
    Int_t mu_charge = 0, ele_charge = 0;
    TauOut mutauh_tau, etauh_tau, tauhtauh_tau1, tauhtauh_tau2;

    // Event bookkeeping and common reconstructed quantities.
    out->Branch("run", &out_run);
    out->Branch("luminosityBlock", &out_lumi);
    out->Branch("event", &out_event);
    out->Branch("genWeight", &gen_weight);
    out->Branch("MET_pt", &met_pt);
    out->Branch("MET_phi", &met_phi);
    out->Branch("HT", &ht);
    out->Branch("n_loose_electron", &n_loose_e);
    out->Branch("n_loose_muon", &n_loose_mu);
    out->Branch("n_boostedTau", &n_boosted_tau);
    out->Branch("n_cleaned_AK4Jet", &n_clean_jet);
    out->Branch("n_cleaned_AK4Jet_btagUParTAK4B_extight", &n_b_extight);
    out->Branch("n_cleaned_AK4Jet_btagUParTAK4B_tight", &n_b_tight);
    out->Branch("n_cleaned_AK4Jet_btagUParTAK4B_medium", &n_b_medium);

    // Generator-level H->tautau labels and four-vectors shared by all channels.
    out->Branch("has_lhe_higgs", &has_lhe_h);
    out->Branch("has_gen_higgs", &has_gen_h);
    out->Branch("is_sig", &is_sig);
    out->Branch("gen_is_mutauh", &is_mutauh);
    out->Branch("gen_is_etauh", &is_etauh);
    out->Branch("gen_is_tauhtauh", &is_tauhtauh);
    out->Branch("gen_decay", &gen_decay);
    setP4Branches(out, "gen_higgs", gen_h_pt, gen_h_eta, gen_h_phi, gen_h_mass);
    setP4Branches(out, "gen_tau1", gen_tau1_pt, gen_tau1_eta, gen_tau1_phi, gen_tau1_mass);
    setP4Branches(out, "gen_tau2", gen_tau2_pt, gen_tau2_eta, gen_tau2_phi, gen_tau2_mass);
    out->Branch("gen_etauh_dr", &gen_etauh_dr);
    out->Branch("gen_mutauh_dr", &gen_mutauh_dr);
    out->Branch("gen_tauhtauh_dr", &gen_tauhtauh_dr);
    out->Branch("gen_tau1_decay", &gen_tau1_decay);
    out->Branch("gen_tau2_decay", &gen_tau2_decay);
    out->Branch("gen_tau1_vis_decayMode", &gen_tau1_vis_dm);
    out->Branch("gen_tau2_vis_decayMode", &gen_tau2_vis_dm);

    // mu tau_h channel: selected muon, selected boostedTau, and matching flags.
    out->Branch("mutauh_reco", &mutauh_reco);
    out->Branch("mutauh_gen_match", &mutauh_gen_match);
    out->Branch("mutauh_higgs_pt", &mutauh_higgs_pt);
    out->Branch("mutauh_dr", &mutauh_dr);
    out->Branch("mutauh_mt", &mutauh_mt);
    setP4Branches(out, "mutauh_muon", mu_pt, mu_eta, mu_phi, mu_mass);
    out->Branch("mutauh_muon_charge", &mu_charge);
    branchTau(out, "mutauh_boostedTau", mutauh_tau);

    // e tau_h channel: selected electron, selected boostedTau, and matching flags.
    out->Branch("etauh_reco", &etauh_reco);
    out->Branch("etauh_gen_match", &etauh_gen_match);
    out->Branch("etauh_higgs_pt", &etauh_higgs_pt);
    out->Branch("etauh_dr", &etauh_dr);
    out->Branch("etauh_mt", &etauh_mt);
    setP4Branches(out, "etauh_electron", ele_pt, ele_eta, ele_phi, ele_mass);
    out->Branch("etauh_electron_charge", &ele_charge);
    branchTau(out, "etauh_boostedTau", etauh_tau);

    // tau_h tau_h channel: selected boostedTau pair and matching flags.
    out->Branch("tauhtauh_reco", &tauhtauh_reco);
    out->Branch("tauhtauh_gen_match", &tauhtauh_gen_match);
    out->Branch("tauhtauh_higgs_pt", &tauhtauh_higgs_pt);
    out->Branch("tauhtauh_dr", &tauhtauh_dr);
    branchTau(out, "tauhtauh_boostedTau1", tauhtauh_tau1);
    branchTau(out, "tauhtauh_boostedTau2", tauhtauh_tau2);

    // Main event loop. Events without a baseline boostedTau object are skipped.
    Long64_t nRead = 0, nFill = 0;
    while (reader.Next()) {
        ++nRead;
        if (nRead % 10000 == 0) std::cout << "Processed " << nRead << " events, filled " << nFill << std::endl;

        // Build the boostedTau baseline collection used by all three channels.
        std::vector<RecoObj> taus;
        for (int i = 0; i < *bt.n; ++i) {
            if (bt.pt[i] > 30.f && std::abs(bt.eta[i]) < 2.4f) {
                taus.push_back({i, bt.pt[i], bt.eta[i], bt.phi[i], bt.mass[i], bt.charge[i]});
            }
        }
        if (taus.empty()) continue;

        // Build common reco and truth objects once per surviving event.
        const auto eles = looseElectrons(lep);
        const auto mus = looseMuons(lep);
        const CommonReco reco = buildCommonReco(jets, eles, mus);
        const GenHTauTau gh = buildGenHTauTau(gen);
        TLorentzVector met;
        met.SetPtEtaPhiM(*ev.metPt, 0.f, *ev.metPhi, 0.f);

        // Reset per-channel outputs so failed channels keep explicit missing values.
        mutauh_reco = mutauh_gen_match = etauh_reco = etauh_gen_match = tauhtauh_reco = tauhtauh_gen_match = 0;
        mutauh_higgs_pt = etauh_higgs_pt = tauhtauh_higgs_pt = kMissing;
        mutauh_dr = etauh_dr = tauhtauh_dr = kMissing;
        mutauh_mt = etauh_mt = kMissing;
        mu_pt = mu_eta = mu_phi = mu_mass = ele_pt = ele_eta = ele_phi = ele_mass = kMissing;
        mu_charge = ele_charge = 0;
        mutauh_tau.reset(); etauh_tau.reset(); tauhtauh_tau1.reset(); tauhtauh_tau2.reset();

        // mu tau_h channel: veto loose electrons and select the best candidate by the configured rule.
        if (eles.empty() && !mus.empty()) {
            int bestMu = -1, bestTau = -1;
            bool bestGenMatched = false;
            for (const auto& mu : mus) {
                for (const auto& tau : taus) {
                    const float dr = deltaR(mu.p4(), tau.p4());
                    if (dr <= 0.1f || dr >= 0.8f || mu.charge * tau.charge >= 0) continue;
                    const float mt = transverseMass(mu.pt, mu.phi, *ev.metPt, *ev.metPhi);
                    if (mt >= 80.f) continue;
                    const float hpt = (mu.p4() + tau.p4() + met).Pt();
                    if (hpt <= 200.f) continue;
                    const bool candGenMatched = gh.isSignal && gh.isMuTauH && gh.hasGenMu &&
                                                deltaR(mu.p4(), gh.genMu.p4()) < 0.4f &&
                                                tauMatchesAnyVisible(tau, gh);
                    if (!preferCandidate(selectionRule, hpt, candGenMatched, mutauh_higgs_pt, bestGenMatched)) continue;
                    mutauh_reco = 1; bestMu = mu.idx; bestTau = tau.idx;
                    mutauh_higgs_pt = hpt; mutauh_dr = dr; mutauh_mt = mt;
                    bestGenMatched = candGenMatched;
                }
            }
            if (mutauh_reco) {
                // Store the best pair and perform gen matching only for true mu tau_h signal events.
                RecoObj mu{bestMu, lep.mPt[bestMu], lep.mEta[bestMu], lep.mPhi[bestMu], lep.mMass[bestMu], lep.mCharge[bestMu]};
                RecoObj tau{bestTau, bt.pt[bestTau], bt.eta[bestTau], bt.phi[bestTau], bt.mass[bestTau], bt.charge[bestTau]};
                fillP4(mu, mu_pt, mu_eta, mu_phi, mu_mass);
                mu_charge = mu.charge;
                fillTau(bt, bestTau, mutauh_tau);
                mutauh_gen_match = bestGenMatched ? 1 : 0;
            }
        }

        // e tau_h channel: veto loose muons and select the best candidate by the configured rule.
        if (mus.empty() && !eles.empty()) {
            int bestEle = -1, bestTau = -1;
            bool bestGenMatched = false;
            for (const auto& ele : eles) {
                for (const auto& tau : taus) {
                    const float dr = deltaR(ele.p4(), tau.p4());
                    if (dr <= 0.1f || dr >= 0.8f || ele.charge * tau.charge >= 0) continue;
                    const float mt = transverseMass(ele.pt, ele.phi, *ev.metPt, *ev.metPhi);
                    if (mt >= 80.f) continue;
                    const float hpt = (ele.p4() + tau.p4() + met).Pt();
                    if (hpt <= 200.f) continue;
                    const bool candGenMatched = gh.isSignal && gh.isETauH && gh.hasGenEle &&
                                                deltaR(ele.p4(), gh.genEle.p4()) < 0.4f &&
                                                tauMatchesAnyVisible(tau, gh);
                    if (!preferCandidate(selectionRule, hpt, candGenMatched, etauh_higgs_pt, bestGenMatched)) continue;
                    etauh_reco = 1; bestEle = ele.idx; bestTau = tau.idx;
                    etauh_higgs_pt = hpt; etauh_dr = dr; etauh_mt = mt;
                    bestGenMatched = candGenMatched;
                }
            }
            if (etauh_reco) {
                // Store the best pair and perform gen matching only for true e tau_h signal events.
                RecoObj ele{bestEle, lep.ePt[bestEle], lep.eEta[bestEle], lep.ePhi[bestEle], lep.eMass[bestEle], lep.eCharge[bestEle]};
                RecoObj tau{bestTau, bt.pt[bestTau], bt.eta[bestTau], bt.phi[bestTau], bt.mass[bestTau], bt.charge[bestTau]};
                fillP4(ele, ele_pt, ele_eta, ele_phi, ele_mass);
                ele_charge = ele.charge;
                fillTau(bt, bestTau, etauh_tau);
                etauh_gen_match = bestGenMatched ? 1 : 0;
            }
        }

        // tau_h tau_h channel: veto loose leptons and select the best opposite-sign pair by the configured rule.
        if (eles.empty() && mus.empty() && taus.size() >= 2) {
            int bestTau1 = -1, bestTau2 = -1;
            bool bestGenMatched = false;
            for (size_t i = 0; i < taus.size(); ++i) {
                for (size_t j = i + 1; j < taus.size(); ++j) {
                    const float dr = deltaR(taus[i].p4(), taus[j].p4());
                    if (dr <= 0.1f || dr >= 0.8f || taus[i].charge * taus[j].charge >= 0) continue;
                    const float hpt = (taus[i].p4() + taus[j].p4() + met).Pt();
                    if (hpt <= 200.f) continue;
                    const bool match12 = gh.hasVisTau1 && gh.hasVisTau2 &&
                                         deltaR(taus[i].p4(), gh.visTau1.p4()) < 0.4f &&
                                         deltaR(taus[j].p4(), gh.visTau2.p4()) < 0.4f;
                    const bool match21 = gh.hasVisTau1 && gh.hasVisTau2 &&
                                         deltaR(taus[i].p4(), gh.visTau2.p4()) < 0.4f &&
                                         deltaR(taus[j].p4(), gh.visTau1.p4()) < 0.4f;
                    const bool candGenMatched = gh.isSignal && gh.isTauHTauH && (match12 || match21);
                    if (!preferCandidate(selectionRule, hpt, candGenMatched, tauhtauh_higgs_pt, bestGenMatched)) continue;
                    tauhtauh_reco = 1; bestTau1 = taus[i].idx; bestTau2 = taus[j].idx;
                    tauhtauh_higgs_pt = hpt; tauhtauh_dr = dr;
                    bestGenMatched = candGenMatched;
                }
            }
            if (tauhtauh_reco) {
                // Store the best pair and its gen matching decision.
                fillTau(bt, bestTau1, tauhtauh_tau1);
                fillTau(bt, bestTau2, tauhtauh_tau2);
                tauhtauh_gen_match = bestGenMatched ? 1 : 0;
            }
        }

        // Fill common event identifiers, reco counts, and truth flags after all channel decisions.
        out_run = *ev.run; out_lumi = *ev.luminosityBlock; out_event = *ev.event;
        gen_weight = *ev.genWeight;
        met_pt = *ev.metPt; met_phi = *ev.metPhi;
        ht = reco.ht; n_loose_e = reco.nLooseEle; n_loose_mu = reco.nLooseMu; n_clean_jet = reco.nCleanJet;
        n_b_extight = reco.nBExTight; n_b_tight = reco.nBTight; n_b_medium = reco.nBMedium;
        n_boosted_tau = taus.size();
        has_lhe_h = gh.hasLHEHiggs; has_gen_h = gh.hasHiggs; is_sig = gh.isSignal;
        is_mutauh = gh.isMuTauH; is_etauh = gh.isETauH; is_tauhtauh = gh.isTauHTauH; gen_decay = gh.decayCode;
        gen_h_pt = gen_h_eta = gen_h_phi = gen_h_mass = kMissing;
        gen_tau1_pt = gen_tau1_eta = gen_tau1_phi = gen_tau1_mass = kMissing;
        gen_tau2_pt = gen_tau2_eta = gen_tau2_phi = gen_tau2_mass = kMissing;
        if (gh.hasHiggs) fillP4(gh.higgs, gen_h_pt, gen_h_eta, gen_h_phi, gen_h_mass);
        if (gh.isSignal) {
            fillP4(gh.tau1, gen_tau1_pt, gen_tau1_eta, gen_tau1_phi, gen_tau1_mass);
            fillP4(gh.tau2, gen_tau2_pt, gen_tau2_eta, gen_tau2_phi, gen_tau2_mass);
        }
        const GenChannelDR ghDR = genChannelDR(gh);
        gen_etauh_dr = ghDR.etauh;
        gen_mutauh_dr = ghDR.mutauh;
        gen_tauhtauh_dr = ghDR.tauhtauh;
        gen_tau1_decay = gh.tau1Decay; gen_tau2_decay = gh.tau2Decay;
        gen_tau1_vis_dm = gh.tau1VisDecayMode; gen_tau2_vis_dm = gh.tau2VisDecayMode;

        // Store one row per event that has at least one baseline boostedTau object.
        out->Fill();
        ++nFill;
    }

    // Persist the tree and close files explicitly for batch jobs.
    writeMetadataTree(inFile, outFile);
    outFile->Write();
    outFile->Close();
    inFile->Close();
    std::cout << "Done. Read " << nRead << " events, filled " << nFill << " events." << std::endl;
}
