//
// Adopted from Erik - https://gitlab.au.dk/ausa/erik/e21010
// Modified with addition of pad spectras


#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1D.h>

using namespace EUtil;
using namespace std;

map<string, vector<string>> det_run_map = { //dssd runs taken before experiment
    {"U1", {"008_000.root", "008_001.root", "008_002.root", "008_003.root"}},
    {"U2", {"008_000.root", "008_001.root", "008_002.root", "008_003.root"}},
    {"U3", {"009_000.root", "009_001.root", "009_002.root", "009_003.root"}},
    {"U4", {"009_000.root", "009_001.root", "009_002.root", "009_003.root"}},
    {"U5", {"010_000.root", "010_001.root", "010_002.root", "010_003.root", "010_004.root", "010_005.root"}},
    {"U6", {"011_000.root", "011_001.root", "011_002.root", "011_003.root", "011_004.root"}}
};/*
map<string, vector<string>> det_run_map = { //dssd runs taken after experiment
    {"U1", {"128_000.root", "128_001.root"}},
    {"U2", {"128_000.root", "128_001.root"}},
    {"U3", {"126_000.root", "126_001.root"}},
    {"U4", {"126_000.root", "126_001.root"}},
    {"U5", {"127_000.root", "127_001.root"}},
    {"U6", {"124_000.root", "124_001.root"}}
};*/
map<string, vector<string>> pad_run_map = { //pad runs taken after experiment
    {"P1", {"129_000.root"}},
    {"P2", {"129_000.root"}},
    {"P3", {"130_000.root"}},
    {"P4", {"130_000.root"}},
    {"P6", {"131_000.root"}}
};
/*
map<string, vector<string>> pad_run_map = { //pad runs taken before experiment
    {"P1", {"005_000.root", "005_001.root"}},
    {"P2", {"003_000.root"}},
    {"P3", {"002_000.root", "002_001.root"}},
    {"P4", {"002_000.root", "002_001.root"}},
    {"P6", {"001_000.root", "001_001.root"}}
};*/


typedef struct detector_side {
  UInt_t *mul_branch{}, *pad_branch{}, *index_branch{}, *energy_branch{};
  string name;
} detector_side;

string input_dir  = getProjectRoot() + "/data/unpacked/calfiles/";
string output_dir = getProjectRoot() + "/data/cal";

int main(int argc, char* argv[]) {
  clock_t start = clock();

  system(("mkdir -p " + output_dir).c_str());
  std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/3ah.root").c_str(), "recreate");

  for (const auto& [det, runs] : det_run_map) {
    auto *c = new TChain("h101");
    for (const auto& run : runs) {
      cout << "Chaining " << run << endl;
      c->Add((input_dir + "/" + run).c_str());
    }

    int strips = 16;
    map<string, TH1D*> hists;
    for (const string& side : {"F", "B"}) {
      for (int i = 0; i < strips; i++) {
        string name = det + side + to_string(i+1);
        Double_t xlow = 1000.; Double_t xup = 2100.;
        if (name == "U6B7") {
          xlow = 500.; xup = 1000.;
        }
        hists.emplace(name, new TH1D(name.c_str(), name.c_str(), (Int_t) (xup - xlow), xlow, xup));
      }
    }

    UInt_t UxF, UxB;
    UInt_t UxF_E[strips], UxFI[strips], UxB_E[strips], UxBI[strips];
    c->SetBranchAddress((det + "F").c_str(), &UxF); c->SetBranchAddress((det + "FI").c_str(), UxFI);
    c->SetBranchAddress((det + "F_E").c_str(), UxF_E);
    c->SetBranchAddress((det + "B").c_str(), &UxB); c->SetBranchAddress((det + "BI").c_str(), UxBI);
    c->SetBranchAddress((det + "B_E").c_str(), UxB_E);

    UInt_t mul, index, energy;
    string name;
    cout << "Total events (" << det << "): " << c->GetEntries() << endl;
    for (int i = 0; i < c->GetEntries(); i++) {
      c->GetEntry(i);

      mul = UxF;
      for (UInt_t j = 0; j < mul; j++) {
        index = UxFI[j];
        energy = UxF_E[j];
        hists.at(det + "F" + to_string(index))->Fill(energy);
      }
      mul = UxB;
      for (UInt_t j = 0; j < mul; j++) {
        index = UxBI[j];
        energy = UxB_E[j];
        hists.at(det + "B" + to_string(index))->Fill(energy);
      }

      if (i % 100000 == 0) cout << "\r" << ceil((double) 100*i/c->GetEntries()) << "% of events processed" << flush;
    }
    cout << endl << endl;

    for (const auto& hist : hists) hist.second->Write();
  }
    // Loop for pad detectors

    cout << "--- Now processing Pad Detectors ---" << endl;
    for (const auto& [pad_name, runs] : pad_run_map) {
        auto* c = new TChain("h101");
        for (const auto& run : runs) {
            cout << "Chaining " << run << " for pad detector " << pad_name << endl;
            c->Add((input_dir + "/" + run).c_str());
        }

        Double_t xlow = 0.; Double_t xup = 4000.; 
        auto* pad_hist = new TH1D(pad_name.c_str(), pad_name.c_str(), (Int_t)xup, xlow, xup);

        UInt_t pad_energy = 0;
        c->SetBranchAddress((pad_name+"E").c_str(), &pad_energy);

        cout << "Total events (" << pad_name << "): " << c->GetEntries() << endl;
        for (int i = 0; i < c->GetEntries(); i++) {
            c->GetEntry(i);

            if (pad_energy > 0) {
                pad_hist->Fill(pad_energy);
            }

             if (i % 100000 == 0) cout << "\r" << ceil((double) 100 * i / c->GetEntries()) << "% of events processed" << flush;
        }
        cout << endl << endl;

        pad_hist->Write();
        delete c;
    }

    out->Close();

    clock_t stop = clock();
    double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
    printf("\nTime elapsed: %.5f seconds\n", elapsed);

    return 0;
}