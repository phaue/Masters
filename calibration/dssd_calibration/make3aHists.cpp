//
// Created by erik on 11/12/23.
//

#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1D.h>

using namespace EUtil;
using namespace std;

map<string, vector<string>> det_run_map = {
    {"U1", {"008_000.root", "008_001.root", "008_002.root", "008_003.root"}},
    {"U2", {"008_000.root", "008_001.root", "008_002.root", "008_003.root"}},
    {"U3", {"009_000.root", "009_001.root", "009_002.root", "009_003.root"}},
    {"U4", {"009_000.root", "009_001.root", "009_002.root", "009_003.root"}},
    {"U5", {"010_000.root", "010_001.root", "010_002.root", "010_003.root", "010_004.root", "010_005.root"}},
    {"U6", {"011_000.root", "011_001.root", "011_002.root", "011_003.root", "011_004.root"}}
};

// TODO: only looking at pre-experiment 3a for now. look also at post-experiment 3a and compare.

typedef struct detector_side {
  UInt_t *mul_branch{}, *pad_branch{}, *index_branch{}, *energy_branch{};
  string name;
} detector_side;

string input_dir  = getProjectRoot() + "/data/unpacked";
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

  out->Close();

  clock_t stop = clock();
  double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
  printf("\nTime elapsed: %.5f seconds\n", elapsed);

  return 0;
}