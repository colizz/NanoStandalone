#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

// Standalone ROOT/ACLiC macro. No generator decay filtering.
void nanoAOD_to_cmsdas_bbvv(const char *input, const char *output = "ntuple.root",
                            double btagWP = 0.1272, double jetPtMin = 15., double jetEtaMax = 2.5,
                            double lepPtMin = 10., double electronIsoMax = 0.12,
                            double muonIsoMax = 0.25) {
    try {
        std::unique_ptr<TFile> fin(TFile::Open(input));
        if (!fin || fin->IsZombie())
            throw std::runtime_error("Cannot open input");
        auto *events = dynamic_cast<TTree *>(fin->Get("Events"));
        if (!events)
            throw std::runtime_error("Missing Events tree");

        TTreeReader r(events);

#define ARRAY(TYPE, NAME) TTreeReaderArray<TYPE> NAME(r, #NAME)
#define VALUE(TYPE, NAME) TTreeReaderValue<TYPE> NAME(r, #NAME)
        ARRAY(Float_t, Jet_pt);
        ARRAY(Float_t, Jet_eta);
        ARRAY(Float_t, Jet_phi);
        ARRAY(Float_t, Jet_mass);
        ARRAY(Float_t, Jet_btagUParTAK4B);

        ARRAY(Float_t, Electron_pt);
        ARRAY(Float_t, Electron_eta);
        ARRAY(Float_t, Electron_phi);
        ARRAY(Float_t, Electron_mass);
        ARRAY(Float_t, Electron_pfRelIso03_all);
        ARRAY(Int_t, Electron_charge);

        ARRAY(Float_t, Muon_pt);
        ARRAY(Float_t, Muon_eta);
        ARRAY(Float_t, Muon_phi);
        ARRAY(Float_t, Muon_mass);
        ARRAY(Float_t, Muon_pfRelIso04_all);
        ARRAY(Int_t, Muon_charge);

        VALUE(Float_t, PuppiMET_pt);
        VALUE(Float_t, PuppiMET_phi);
        VALUE(Float_t, genWeight);
        VALUE(UInt_t, run);
        VALUE(UInt_t, luminosityBlock);
        VALUE(ULong64_t, event);
#undef ARRAY
#undef VALUE

        // Fail even for empty inputs if a required branch is absent.
        for (const char *n : {"Jet_pt",
                              "Jet_eta",
                              "Jet_phi",
                              "Jet_mass",
                              "Jet_btagUParTAK4B",
                              "Electron_pt",
                              "Electron_eta",
                              "Electron_phi",
                              "Electron_mass",
                              "Electron_pfRelIso03_all",
                              "Electron_charge",
                              "Muon_pt",
                              "Muon_eta",
                              "Muon_phi",
                              "Muon_mass",
                              "Muon_pfRelIso04_all",
                              "Muon_charge",
                              "PuppiMET_pt",
                              "PuppiMET_phi",
                              "genWeight",
                              "run",
                              "luminosityBlock",
                              "event"})
            if (!events->GetBranch(n))
                throw std::runtime_error(std::string("Missing branch: ") + n);

        TFile fout(output, "RECREATE");
        if (fout.IsZombie())
            throw std::runtime_error("Cannot create output");

        TTree tree("tree", "CMSDAS compatible objects");
#define VF(N)                                                                                      \
    std::vector<float> N;                                                                          \
    tree.Branch(#N, &N)
#define VI(N)                                                                                      \
    std::vector<int> N;                                                                            \
    tree.Branch(#N, &N)
        VF(jet_pt);
        VF(jet_eta);
        VF(jet_phi);
        VF(jet_energy);
        VF(jet_mass);
        VI(jet_is_btag);
        VF(lep_pt);
        VF(lep_eta);
        VF(lep_phi);
        VF(lep_energy);
        VI(lep_charge);
        VI(lep_pid);
        VF(lep_iso);
#undef VF
#undef VI

        float met_pt = 0, met_phi = 0;
        tree.Branch("met_pt", &met_pt, "met_pt/F");
        tree.Branch("met_phi", &met_phi, "met_phi/F");

        TTree audit("audit", "One row per selected event; aligned with tree");
        UInt_t outRun = 0, outLumi = 0;
        ULong64_t outEvent = 0;
        float weight = 0;
        audit.Branch("run", &outRun, "run/i");
        audit.Branch("luminosityBlock", &outLumi, "luminosityBlock/i");
        audit.Branch("event", &outEvent, "event/l");
        audit.Branch("genWeight", &weight, "genWeight/F");

        TH1D cutflow("cutflow", "Unweighted event counts", 3, 0, 3);
        TH1D weights("sumGenWeight", "Generator weight sums at each stage", 3, 0, 3);
        weights.Sumw2();
        const char *labels[] = {"input", "at_least_2_btags", "at_least_2_leptons"};
        for (int i = 0; i < 3; ++i) {
            cutflow.GetXaxis()->SetBinLabel(i + 1, labels[i]);
            weights.GetXaxis()->SetBinLabel(i + 1, labels[i]);
        }

        auto energy = [](double pt, double eta, double m) {
            return static_cast<float>(std::hypot(pt * std::cosh(eta), m));
        };

        Long64_t processed = 0;

        while (r.Next()) {
            ++processed;
            cutflow.Fill(0.5);
            weights.Fill(0.5, *genWeight);

            jet_pt.clear();
            jet_eta.clear();
            jet_phi.clear();
            jet_energy.clear();
            jet_mass.clear();
            jet_is_btag.clear();

            lep_pt.clear();
            lep_eta.clear();
            lep_phi.clear();
            lep_energy.clear();
            lep_charge.clear();
            lep_pid.clear();
            lep_iso.clear();

            int nb = 0;
            for (size_t i = 0; i < Jet_pt.GetSize(); ++i) {
                if (Jet_pt[i] <= jetPtMin || std::abs(Jet_eta[i]) >= jetEtaMax)
                    continue;
                jet_pt.push_back(Jet_pt[i]);
                jet_eta.push_back(Jet_eta[i]);
                jet_phi.push_back(Jet_phi[i]);
                jet_mass.push_back(Jet_mass[i]);
                jet_energy.push_back(energy(Jet_pt[i], Jet_eta[i], Jet_mass[i]));
                int tagged = Jet_btagUParTAK4B[i] > btagWP;
                jet_is_btag.push_back(tagged);
                nb += tagged;
            }
            if (nb < 2)
                continue;
            cutflow.Fill(1.5);
            weights.Fill(1.5, *genWeight);

            auto add =
                [&](float pt, float eta, float phi, float mass, int charge, int absPid, float iso) {
                    lep_pt.push_back(pt);
                    lep_eta.push_back(eta);
                    lep_phi.push_back(phi);
                    lep_energy.push_back(energy(pt, eta, mass));
                    lep_charge.push_back(charge);
                    lep_pid.push_back(-absPid * charge);
                    lep_iso.push_back(iso);
                };

            for (size_t i = 0; i < Electron_pt.GetSize(); ++i)
                if (Electron_pt[i] > lepPtMin && std::abs(Electron_eta[i]) < 2.5 &&
                    Electron_pfRelIso03_all[i] < electronIsoMax)
                    add(Electron_pt[i],
                        Electron_eta[i],
                        Electron_phi[i],
                        Electron_mass[i],
                        Electron_charge[i],
                        11,
                        Electron_pfRelIso03_all[i]);

            for (size_t i = 0; i < Muon_pt.GetSize(); ++i)
                if (Muon_pt[i] > lepPtMin && std::abs(Muon_eta[i]) < 2.4 &&
                    Muon_pfRelIso04_all[i] < muonIsoMax)
                    add(Muon_pt[i],
                        Muon_eta[i],
                        Muon_phi[i],
                        Muon_mass[i],
                        Muon_charge[i],
                        13,
                        Muon_pfRelIso04_all[i]);
            if (lep_pt.size() < 2)
                continue;
            cutflow.Fill(2.5);
            weights.Fill(2.5, *genWeight);
            met_pt = *PuppiMET_pt;
            met_phi = *PuppiMET_phi;
            outRun = *run;
            outLumi = *luminosityBlock;
            outEvent = *event;
            weight = *genWeight;
            if (tree.Fill() < 0 || audit.Fill() < 0)
                throw std::runtime_error("Output Fill failed");
        }

        if (processed != events->GetEntries())
            throw std::runtime_error("Incomplete read (branch type or I/O error)");

        TNamed source("input_source", input);
        source.Write();

        TNamed config("selection",
                      Form("Jet pt>%g abs(eta)<%g; UParTAK4B>%g; leptons pt>%g; electron "
                           "abs(eta)<2.5 iso03<%g; muon abs(eta)<2.4 iso04<%g; >=2 btags and >=2 "
                           "leptons; no ID, overlap, trigger or truth cuts",
                           jetPtMin,
                           jetEtaMax,
                           btagWP,
                           lepPtMin,
                           electronIsoMax,
                           muonIsoMax));
        config.Write();
        if (fout.Write() <= 0 || fout.TestBit(TFile::kWriteError))
            throw std::runtime_error("Output Write failed");
        std::cout << "Processed " << processed << ", selected " << tree.GetEntries() << std::endl;
        fout.Close();
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        gSystem->Exit(1);
    }
}
