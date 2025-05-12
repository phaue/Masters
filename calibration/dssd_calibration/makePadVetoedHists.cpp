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

vector<string> files = {
    "098_000.root", "098_001.root", "098_002.root",
    "099_000.root", "099_001.root", "099_002.root",
    "100_000.root", "100_001.root", "100_002.root",
    "101_000.root", "101_001.root"
};

typedef struct detector_side {
  UInt_t *mul_branch{}, *pad_branch{}, *index_branch{}, *energy_branch{};
  string name;
} detector_side;

string input_dir  = getProjectRoot() + "/data/unpacked/Si";
string output_dir = getProjectRoot() + "/data/cal";

int main(int argc, char* argv[]) {
  clock_t start = clock();

  auto *c = new TChain("h101");
  for (const auto& file : files) {
    cout << "Chaining " << file << endl;
    c->Add((input_dir + "/" + file).c_str());
  }

  int strips = 16;
  Double_t xlow = 50.; Double_t xup = 1000.;
  map<string, TH1D*> hists;
  for (const string& num : {"1", "2", "3"}) { // U4 is 300 um, U5 is 1000 um - no pad vetoing
    for (const string& side : {"F", "B"}) {
      for (int i = 0; i < strips; i++) {
        string name = "U" + num + side + to_string(i+1);
        hists.emplace(name, new TH1D(name.c_str(), name.c_str(), (Int_t) (xup - xlow), xlow, xup));
      }
    }
  }

  UInt_t U1F, U1B, U2F, U2B, U3F, U3B, P1E, P2E, P3E;
  UInt_t U1F_E[strips], U1FI[strips], U1B_E[strips], U1BI[strips];
  UInt_t U2F_E[strips], U2FI[strips], U2B_E[strips], U2BI[strips];
  UInt_t U3F_E[strips], U3FI[strips], U3B_E[strips], U3BI[strips];
  c->SetBranchAddress("U1F", &U1F); c->SetBranchAddress("U1FI", U1FI); c->SetBranchAddress("U1F_E", U1F_E);
  c->SetBranchAddress("U1B", &U1B); c->SetBranchAddress("U1BI", U1BI); c->SetBranchAddress("U1B_E", U1B_E);
  c->SetBranchAddress("U2F", &U2F); c->SetBranchAddress("U2FI", U2FI); c->SetBranchAddress("U2F_E", U2F_E);
  c->SetBranchAddress("U2B", &U2B); c->SetBranchAddress("U2BI", U2BI); c->SetBranchAddress("U2B_E", U2B_E);
  c->SetBranchAddress("U3F", &U3F); c->SetBranchAddress("U3FI", U3FI); c->SetBranchAddress("U3F_E", U3F_E);
  c->SetBranchAddress("U3B", &U3B); c->SetBranchAddress("U3BI", U3BI); c->SetBranchAddress("U3B_E", U3B_E);
  c->SetBranchAddress("P1E", &P1E); c->SetBranchAddress("P2E", &P2E);
  c->SetBranchAddress("P3E", &P3E);

  system(("mkdir -p " + output_dir).c_str());
  std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/pvh.root").c_str(), "recreate");
//  auto *h1 = new TH1D("h1", "h1", 4000, 0, 4000);

  vector<detector_side> detector_sides = {
      {&U1F, &P1E, U1FI, U1F_E, "U1F"},
      {&U1B, &P1E, U1BI, U1B_E, "U1B"},
      {&U2F, &P2E, U2FI, U2F_E, "U2F"},
      {&U2B, &P2E, U2BI, U2B_E, "U2B"},
      {&U3F, &P3E, U3FI, U3F_E, "U3F"},
      {&U3B, &P3E, U3BI, U3B_E, "U3B"},
  };
  UInt_t mul, pad_energy, index, energy;
  string name;
  cout << "Total events: " << c->GetEntries() << endl;
  for (int i = 0; i < c->GetEntries(); i++) {
    c->GetEntry(i);

    for (const auto& side : detector_sides) {
      pad_energy = *side.pad_branch;
      if (pad_energy == 0) {
        mul = *side.mul_branch;
        name = side.name;
        for (UInt_t j = 0; j < mul; j++) {
          index = side.index_branch[j];
          energy = side.energy_branch[j];
          hists.at(name + to_string(index))->Fill(energy);
        }
      }
    }
    if (i % 100000 == 0) cout << "\r" << ceil((double) 100*i/c->GetEntries()) << "% of events processed" << flush;
  }
  cout << endl;

  for (const auto& hist : hists) hist.second->Write();

  out->Close();

  clock_t stop = clock();
  double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
  printf("\nTime elapsed: %.5f seconds\n", elapsed);

  return 0;
}