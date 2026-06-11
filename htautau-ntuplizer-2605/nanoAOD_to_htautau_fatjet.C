/*
 * AK8 fatjet H->tautau ntuplizer.
 *
 * Purpose:
 *   Build a compact event-level ntuple for the fatjet-based boosted H->tautau
 *   reconstruction strategy. Each selected event contributes at most one AK8
 *   candidate, plus common event, reco-summary, generator-truth, and metadata
 *   information. The output can be compared directly with the boostedTau
 *   ntuplizer output.
 *
 * Input and output:
 *   - Reads NanoAOD Events and Runs trees. Plain /store paths are resolved
 *     through the XRootD redirector in the shared utility header.
 *   - Writes an Events-like tree named HTauTauFatJet.
 *   - Writes a Metadata tree from Runs with genEventCount, total gen weight
 *     / genEventSumw, and genEventSumw2 when available.
 *   - Per-event genWeight from Events is stored in the main ntuple.
 *
 * Event and object selections:
 *   - Events with no NanoAOD FatJet objects are skipped immediately.
 *   - A baseline AK8 candidate must satisfy raw pT = pt*(1-rawFactor) > 200 GeV
 *     and |eta| < 2.4. Events without any such baseline AK8 jet are skipped.
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
 * Fatjet candidate rule:
 *   object_selection_rule is exposed as the third macro argument.
 *   - "max_pt" (default): choose the baseline AK8 jet with largest raw pT.
 *   - "max_tau_score": choose the baseline AK8 jet with largest
 *     (GlobalParT3_Xtauhtaue + Xtauhtaum + Xtauhtauh)
 *     / (same tau score + GlobalParT3_QCD).
 *   - "gen_matched": choose the baseline AK8 jet with dR(gen Higgs, jet) < 0.8
 *     and smallest dR. If no such jet exists, the event is still kept when it
 *     has a baseline AK8 jet, but selected-fatjet variables are filled with
 *     dummy values.
 *
 * Generator truth and matching:
 *   - Finds the last-copy gen Higgs from GenPart and labels H->tautau signal.
 *   - Tau decays are classified as e tau_h, mu tau_h, tau_h tau_h, or other.
 *   - LHEPart is optional: has_lhe_higgs is filled only when LHE branches exist.
 *   - gen_etauh_dr, gen_mutauh_dr, and gen_tauhtauh_dr are filled only for the
 *     corresponding gen channel using gen e/mu and visible hadronic tau objects;
 *     otherwise they are set to kMissing.
 *   - fatjet_gen_match is true when the selected AK8 jet is within dR < 0.8
 *     of the gen Higgs.
 *
 * Stored variable groups:
 *   - Event: run, luminosityBlock, event, genWeight, MET_pt, MET_phi.
 *   - Reco summary: HT, loose lepton counts, boostedTau count, cleaned AK4 jet
 *     count, cleaned AK4 b-tag counts at medium/tight/ex-tight WPs.
 *   - Gen: Higgs/tau four-vectors, decay labels, tau decay modes, visible tau
 *     decay modes, channel dR values, and signal/channel flags.
 *   - FatJet: selected raw pt/mass/soft-drop mass, eta, phi, dR to gen Higgs,
 *     dR to nearest loose electron/muon, GlobalParT3 scores and mass
 *     corrections, ParticleNet scores and mass corrections.
 *
 * Notes:
 *   - Stored AK8 pt, mass, and soft-drop mass are raw-corrected with
 *     (1 - FatJet_rawFactor), while eta and phi are unchanged.
 *   - Dummy float values use kMissing from HTauTauNtuplizerUtils.h.
 */
#include "HTauTauNtuplizerUtils.h"

#include "TFile.h"

// Reader bundle for AK8 kinematics and the H->tautau-oriented tagger scores.
struct FatJetBranches {
    TTreeReaderValue<Int_t> nFatJet;
    TTreeReaderArray<Float_t> pt, eta, phi, mass, msoftdrop, rawFactor;
    TTreeReaderArray<Float_t> globalParT3_QCD, globalParT3_TopbWev, globalParT3_TopbWmv;
    TTreeReaderArray<Float_t> globalParT3_TopbWq, globalParT3_TopbWqq, globalParT3_TopbWtauhv;
    TTreeReaderArray<Float_t> globalParT3_Xqq, globalParT3_Xtauhtaue, globalParT3_Xtauhtauh, globalParT3_Xtauhtaum;
    TTreeReaderArray<Float_t> globalParT3_massCorrGeneric, globalParT3_massCorrX2p;
    TTreeReaderArray<Float_t> particleNetLegacy_mass, particleNet_QCD, particleNet_XqqVsQCD;
    TTreeReaderArray<Float_t> particleNet_XteVsQCD, particleNet_XtmVsQCD, particleNet_XttVsQCD, particleNet_massCorr;

    // Bind only the FatJet branches that are written to the output ntuple.
    FatJetBranches(TTreeReader& r)
        : nFatJet(r, "nFatJet"), pt(r, "FatJet_pt"), eta(r, "FatJet_eta"), phi(r, "FatJet_phi"),
          mass(r, "FatJet_mass"), msoftdrop(r, "FatJet_msoftdrop"), rawFactor(r, "FatJet_rawFactor"),
          globalParT3_QCD(r, "FatJet_globalParT3_QCD"),
          globalParT3_TopbWev(r, "FatJet_globalParT3_TopbWev"),
          globalParT3_TopbWmv(r, "FatJet_globalParT3_TopbWmv"),
          globalParT3_TopbWq(r, "FatJet_globalParT3_TopbWq"),
          globalParT3_TopbWqq(r, "FatJet_globalParT3_TopbWqq"),
          globalParT3_TopbWtauhv(r, "FatJet_globalParT3_TopbWtauhv"),
          globalParT3_Xqq(r, "FatJet_globalParT3_Xqq"),
          globalParT3_Xtauhtaue(r, "FatJet_globalParT3_Xtauhtaue"),
          globalParT3_Xtauhtauh(r, "FatJet_globalParT3_Xtauhtauh"),
          globalParT3_Xtauhtaum(r, "FatJet_globalParT3_Xtauhtaum"),
          globalParT3_massCorrGeneric(r, "FatJet_globalParT3_massCorrGeneric"),
          globalParT3_massCorrX2p(r, "FatJet_globalParT3_massCorrX2p"),
          particleNetLegacy_mass(r, "FatJet_particleNetLegacy_mass"),
          particleNet_QCD(r, "FatJet_particleNet_QCD"),
          particleNet_XqqVsQCD(r, "FatJet_particleNet_XqqVsQCD"),
          particleNet_XteVsQCD(r, "FatJet_particleNet_XteVsQCD"),
          particleNet_XtmVsQCD(r, "FatJet_particleNet_XtmVsQCD"),
          particleNet_XttVsQCD(r, "FatJet_particleNet_XttVsQCD"),
          particleNet_massCorr(r, "FatJet_particleNet_massCorr") {}
};

// GlobalParT3 H->tautau-vs-QCD score used by the max_tau_score option.
float fatJetTauVsQCD(const FatJetBranches& fat, int idx) {
    const float tauScore = fat.globalParT3_Xtauhtaue[idx] + fat.globalParT3_Xtauhtaum[idx] + fat.globalParT3_Xtauhtauh[idx];
    const float denom = tauScore + fat.globalParT3_QCD[idx];
    return denom > 0.f ? tauScore / denom : -1.f;
}

// Convert NanoAOD events into an event-level ntuple for the leading AK8 strategy.
void nanoAOD_to_htautau_fatjet(const char* inputFile, const char* outputFile, const char* object_selection_rule = "max_pt") {
    const std::string selectionRule(object_selection_rule);
    if (selectionRule != "max_pt" && selectionRule != "max_tau_score" && selectionRule != "gen_matched") {
        std::cerr << "Unsupported object_selection_rule=\"" << selectionRule
                  << "\". Supported values are max_pt, max_tau_score, and gen_matched." << std::endl;
        return;
    }

    // Resolve /store paths through XRootD, then open the NanoAOD input.
    const std::string inputPath = inputPathWithRedirector(inputFile);
    TFile* inFile = TFile::Open(inputPath.c_str(), "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error opening input file: " << inputFile << std::endl;
        return;
    }

    // Bind all input collections through TTreeReader wrappers.
    TTreeReader reader("Events", inFile);
    EventBranches ev(reader);
    GenBranches gen(reader);
    LeptonBranches lep(reader);
    JetBranches jets(reader);
    FatJetBranches fat(reader);
    TTreeReaderValue<Int_t> nInputBoostedTau(reader, "nboostedTau");

    // Create one output row per selected event; gen_matched can keep a row with dummy AK8 values.
    TFile* outFile = new TFile(outputFile, "RECREATE");
    TTree* out = new TTree("HTauTauFatJet", "H->tautau AK8 fatjet ntuple");

    // Output storage: event IDs, MET, common reco summary, truth labels, and AK8 payload.
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
    Int_t fatjet_gen_match = 0;
    Float_t fatjet_dr_gen_h = 99.f;
    Float_t fatjet_dr_nearest_loose_electron = 99.f, fatjet_dr_nearest_loose_muon = 99.f;
    Float_t fj_pt = kMissing, fj_eta = kMissing, fj_phi = kMissing, fj_mass = kMissing, fj_sdmass = kMissing;
    Float_t gp_qcd = kMissing, gp_top_e = kMissing, gp_top_m = kMissing, gp_top_q = kMissing, gp_top_qq = kMissing, gp_top_tauh = kMissing;
    Float_t gp_xqq = kMissing, gp_xtauhe = kMissing, gp_xtauhh = kMissing, gp_xtauhm = kMissing, gp_mass_corr_generic = kMissing, gp_mass_corr_x2p = kMissing;
    Float_t pn_legacy_mass = kMissing, pn_qcd = kMissing, pn_xqq = kMissing, pn_xte = kMissing, pn_xtm = kMissing, pn_xtt = kMissing, pn_mass_corr = kMissing;

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

    // Generator-level H->tautau labels and four-vectors.
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

    // Leading AK8 kinematics, raw masses, matching flag, and requested tagger scores.
    out->Branch("fatjet_gen_match", &fatjet_gen_match);
    out->Branch("fatjet_dr_gen_higgs", &fatjet_dr_gen_h);
    out->Branch("fatjet_dr_nearest_loose_electron", &fatjet_dr_nearest_loose_electron);
    out->Branch("fatjet_dr_nearest_loose_muon", &fatjet_dr_nearest_loose_muon);
    out->Branch("fatjet_raw_pt", &fj_pt);
    out->Branch("fatjet_eta", &fj_eta);
    out->Branch("fatjet_phi", &fj_phi);
    out->Branch("fatjet_raw_mass", &fj_mass);
    out->Branch("fatjet_raw_sdmass", &fj_sdmass);
    out->Branch("FatJet_globalParT3_QCD", &gp_qcd);
    out->Branch("FatJet_globalParT3_TopbWev", &gp_top_e);
    out->Branch("FatJet_globalParT3_TopbWmv", &gp_top_m);
    out->Branch("FatJet_globalParT3_TopbWq", &gp_top_q);
    out->Branch("FatJet_globalParT3_TopbWqq", &gp_top_qq);
    out->Branch("FatJet_globalParT3_TopbWtauhv", &gp_top_tauh);
    out->Branch("FatJet_globalParT3_Xqq", &gp_xqq);
    out->Branch("FatJet_globalParT3_Xtauhtaue", &gp_xtauhe);
    out->Branch("FatJet_globalParT3_Xtauhtauh", &gp_xtauhh);
    out->Branch("FatJet_globalParT3_Xtauhtaum", &gp_xtauhm);
    out->Branch("FatJet_globalParT3_massCorrGeneric", &gp_mass_corr_generic);
    out->Branch("FatJet_globalParT3_massCorrX2p", &gp_mass_corr_x2p);
    out->Branch("FatJet_particleNetLegacy_mass", &pn_legacy_mass);
    out->Branch("FatJet_particleNet_QCD", &pn_qcd);
    out->Branch("FatJet_particleNet_XqqVsQCD", &pn_xqq);
    out->Branch("FatJet_particleNet_XteVsQCD", &pn_xte);
    out->Branch("FatJet_particleNet_XtmVsQCD", &pn_xtm);
    out->Branch("FatJet_particleNet_XttVsQCD", &pn_xtt);
    out->Branch("FatJet_particleNet_massCorr", &pn_mass_corr);

    // Main event loop. Events without an AK8 candidate are intentionally skipped.
    Long64_t nRead = 0, nFill = 0;
    Long64_t nFilledWithGenHiggs = 0;
    Long64_t nFilledWithGenHiggsNoFatJetDR08 = 0;
    Long64_t nFilledETauH = 0, nFilledMuTauH = 0, nFilledTauHTauH = 0;
    Long64_t nNoFatJetDR08ETauH = 0, nNoFatJetDR08MuTauH = 0, nNoFatJetDR08TauHTauH = 0;
    while (reader.Next()) {
        ++nRead;
        if (nRead % 10000 == 0) std::cout << "Processed " << nRead << " events, filled " << nFill << std::endl;
        if (*fat.nFatJet == 0) continue;

        // Build common loose leptons, cleaned AK4 summary variables, and truth labels.
        const auto eles = looseElectrons(lep);
        const auto mus = looseMuons(lep);
        const CommonReco reco = buildCommonReco(jets, eles, mus);
        const GenHTauTau gh = buildGenHTauTau(gen);

        // Select the AK8 candidate according to the requested object_selection_rule.
        // In gen_matched mode, keep the event with dummy AK8 values when no jet matches gen Higgs.
        int lead = -1;
        float leadRawPt = kMissing;
        float bestMetric = -1.f;
        float bestDR = 99.f;
        bool hasBaselineFatJet = false;
        for (int i = 0; i < *fat.nFatJet; ++i) {
            const float rawPt = fat.pt[i] * (1.f - fat.rawFactor[i]);
            if (rawPt <= 200.f || std::abs(fat.eta[i]) >= 2.4f) continue;
            hasBaselineFatJet = true;

            if (selectionRule == "max_pt") {
                if (rawPt > bestMetric) {
                    lead = i;
                    leadRawPt = rawPt;
                    bestMetric = rawPt;
                }
            } else if (selectionRule == "max_tau_score") {
                const float tauVsQCD = fatJetTauVsQCD(fat, i);
                if (tauVsQCD > bestMetric) {
                    lead = i;
                    leadRawPt = rawPt;
                    bestMetric = tauVsQCD;
                }
            } else if (selectionRule == "gen_matched" && gh.hasHiggs) {
                TLorentzVector fjForDR;
                fjForDR.SetPtEtaPhiM(rawPt, fat.eta[i], fat.phi[i], fat.mass[i] * (1.f - fat.rawFactor[i]));
                const float dr = deltaR(fjForDR, gh.higgs.p4());
                if (dr < 0.8f && dr < bestDR) {
                    lead = i;
                    leadRawPt = rawPt;
                    bestDR = dr;
                }
            }
        }
        if (!hasBaselineFatJet) continue;
        if (selectionRule != "gen_matched" && lead < 0) continue;

        // For stored events with a gen Higgs, count cases with no selected-quality AK8 jet near it.
        if (gh.hasHiggs) {
            bool hasFatJetNearGenHiggs = false;
            for (int i = 0; i < *fat.nFatJet; ++i) {
                const float rawPt = fat.pt[i] * (1.f - fat.rawFactor[i]);
                if (rawPt <= 200.f || std::abs(fat.eta[i]) >= 2.4f) continue;
                TLorentzVector fjForDR;
                fjForDR.SetPtEtaPhiM(rawPt, fat.eta[i], fat.phi[i], fat.mass[i] * (1.f - fat.rawFactor[i]));
                if (deltaR(fjForDR, gh.higgs.p4()) < 0.8f) {
                    hasFatJetNearGenHiggs = true;
                    break;
                }
            }
            ++nFilledWithGenHiggs;
            if (!hasFatJetNearGenHiggs) ++nFilledWithGenHiggsNoFatJetDR08;
            if (gh.isETauH) {
                ++nFilledETauH;
                if (!hasFatJetNearGenHiggs) ++nNoFatJetDR08ETauH;
            }
            if (gh.isMuTauH) {
                ++nFilledMuTauH;
                if (!hasFatJetNearGenHiggs) ++nNoFatJetDR08MuTauH;
            }
            if (gh.isTauHTauH) {
                ++nFilledTauHTauH;
                if (!hasFatJetNearGenHiggs) ++nNoFatJetDR08TauHTauH;
            }
        }

        // Fill event-level identifiers and reco summary branches.
        out_run = *ev.run; out_lumi = *ev.luminosityBlock; out_event = *ev.event;
        gen_weight = *ev.genWeight;
        met_pt = *ev.metPt; met_phi = *ev.metPhi;
        ht = reco.ht; n_loose_e = reco.nLooseEle; n_loose_mu = reco.nLooseMu; n_clean_jet = reco.nCleanJet;
        n_b_extight = reco.nBExTight; n_b_tight = reco.nBTight; n_b_medium = reco.nBMedium;
        n_boosted_tau = *nInputBoostedTau;

        // Fill truth flags and reset missing truth four-vectors by default.
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

        // Reset AK8 payload. gen_matched events with no matched AK8 keep these dummy values.
        fatjet_gen_match = 0;
        fatjet_dr_gen_h = 99.f;
        fatjet_dr_nearest_loose_electron = 99.f;
        fatjet_dr_nearest_loose_muon = 99.f;
        fj_pt = fj_eta = fj_phi = fj_mass = fj_sdmass = kMissing;
        gp_qcd = gp_top_e = gp_top_m = gp_top_q = gp_top_qq = gp_top_tauh = kMissing;
        gp_xqq = gp_xtauhe = gp_xtauhh = gp_xtauhm = gp_mass_corr_generic = gp_mass_corr_x2p = kMissing;
        pn_legacy_mass = pn_qcd = pn_xqq = pn_xte = pn_xtm = pn_xtt = pn_mass_corr = kMissing;

        if (lead >= 0) {
            // Convert AK8 pT/mass/soft-drop mass back to raw values and match to gen Higgs.
            const float rawScale = 1.f - fat.rawFactor[lead];
            fj_pt = leadRawPt;
            fj_eta = fat.eta[lead];
            fj_phi = fat.phi[lead];
            fj_mass = fat.mass[lead] * rawScale;
            fj_sdmass = fat.msoftdrop[lead] * rawScale;
            TLorentzVector fj;
            fj.SetPtEtaPhiM(fj_pt, fj_eta, fj_phi, fj_mass);
            fatjet_dr_gen_h = gh.hasHiggs ? deltaR(fj, gh.higgs.p4()) : 99.f;
            fatjet_gen_match = (gh.hasHiggs && fatjet_dr_gen_h < 0.8f) ? 1 : 0;
            for (const auto& ele : eles) {
                fatjet_dr_nearest_loose_electron = std::min(fatjet_dr_nearest_loose_electron, deltaR(fj, ele.p4()));
            }
            for (const auto& mu : mus) {
                fatjet_dr_nearest_loose_muon = std::min(fatjet_dr_nearest_loose_muon, deltaR(fj, mu.p4()));
            }

            // Copy the selected AK8 GlobalParT3 and ParticleNet scores.
            gp_qcd = fat.globalParT3_QCD[lead];
            gp_top_e = fat.globalParT3_TopbWev[lead];
            gp_top_m = fat.globalParT3_TopbWmv[lead];
            gp_top_q = fat.globalParT3_TopbWq[lead];
            gp_top_qq = fat.globalParT3_TopbWqq[lead];
            gp_top_tauh = fat.globalParT3_TopbWtauhv[lead];
            gp_xqq = fat.globalParT3_Xqq[lead];
            gp_xtauhe = fat.globalParT3_Xtauhtaue[lead];
            gp_xtauhh = fat.globalParT3_Xtauhtauh[lead];
            gp_xtauhm = fat.globalParT3_Xtauhtaum[lead];
            gp_mass_corr_generic = fat.globalParT3_massCorrGeneric[lead];
            gp_mass_corr_x2p = fat.globalParT3_massCorrX2p[lead];
            pn_legacy_mass = fat.particleNetLegacy_mass[lead];
            pn_qcd = fat.particleNet_QCD[lead];
            pn_xqq = fat.particleNet_XqqVsQCD[lead];
            pn_xte = fat.particleNet_XteVsQCD[lead];
            pn_xtm = fat.particleNet_XtmVsQCD[lead];
            pn_xtt = fat.particleNet_XttVsQCD[lead];
            pn_mass_corr = fat.particleNet_massCorr[lead];
        }

        // Store exactly one row for the selected event.
        out->Fill();
        ++nFill;
    }

    // Persist the tree and close files explicitly for batch jobs.
    writeMetadataTree(inFile, outFile);
    outFile->Write();
    outFile->Close();
    inFile->Close();
    std::cout << "Done. Read " << nRead << " events, filled " << nFill << " events." << std::endl;
    std::cout << "Stored events with gen Higgs: " << nFilledWithGenHiggs << std::endl;
    std::cout << "Stored events with no selected-quality AK8 fatjet within dR<0.8 of gen Higgs: "
              << nFilledWithGenHiggsNoFatJetDR08 << std::endl;
    std::cout << "Stored gen e tau_h events: " << nFilledETauH
              << ", no fatjet within dR<0.8: " << nNoFatJetDR08ETauH << std::endl;
    std::cout << "Stored gen mu tau_h events: " << nFilledMuTauH
              << ", no fatjet within dR<0.8: " << nNoFatJetDR08MuTauH << std::endl;
    std::cout << "Stored gen tau_h tau_h events: " << nFilledTauHTauH
              << ", no fatjet within dR<0.8: " << nNoFatJetDR08TauHTauH << std::endl;
}
