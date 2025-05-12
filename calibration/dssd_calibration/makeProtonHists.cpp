//
// Created by erik on 11/12/23.
// changed to fit new calibrations on 30/04 by philip
//
/*Denne fil tager de unpackede filer og danner histogrammer med de ukalibrerede data fra DSSSD 4 
den laver et histogram pr bagstrip og pr forstrip, dvs 32 histogrammer i alt - det er så dem vi kan lave
videre analyse på i næste script
*/
#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1D.h>

using namespace EUtil;
using namespace std;

vector<string> files = {
    "098_000.root", "098_001.root", "098_002.root", "099_000.root", "099_001.root",
    "099_002.root", "100_000.root", "100_001.root", "100_002.root", "101_000.root",
    "101_001.root"
};

typedef struct detector_side {
  UInt_t *mul_branch{}, *index_branch{}, *energy_branch{};
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
  Double_t xlow = 90.; Double_t xup = 2000.;
  map<string, TH1D*> hists;
  for (const string& num : {"4"}) {
    if (num == "4") xlow = 90.;
    for (const string& side : {"F", "B"}) {
      for (int i = 0; i < strips; i++) {
        string name = "U" + num + side + to_string(i+1);
        hists.emplace(name, new TH1D(name.c_str(), name.c_str(), (Int_t) (xup - xlow), xlow, xup));
      }
    }
  }

  UInt_t U4F, U4B;
  UInt_t U4F_E[strips], U4FI[strips], U4B_E[strips], U4BI[strips];
  c->SetBranchAddress("U4F", &U4F); c->SetBranchAddress("U4FI", U4FI); c->SetBranchAddress("U4F_E", U4F_E);
  c->SetBranchAddress("U4B", &U4B); c->SetBranchAddress("U4BI", U4BI); c->SetBranchAddress("U4B_E", U4B_E);

  system(("mkdir -p " + output_dir).c_str());
  std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/ph.root").c_str(), "recreate");
//  auto *h1 = new TH1D("h1", "h1", 4000, 0, 4000);

  vector<detector_side> detector_sides = {
      {&U4F, U4FI, U4F_E, "U4F"},
      {&U4B, U4BI, U4B_E, "U4B"},
  };
  UInt_t mul, index, energy;
  string name;
  cout << "Total events: " << c->GetEntries() << endl;
  for (int i = 0; i < c->GetEntries(); i++) {
    c->GetEntry(i);

    for (const auto& side : detector_sides) {
      mul = *side.mul_branch;
      name = side.name;
      for (UInt_t j = 0; j < mul; j++) {
        index = side.index_branch[j];
        energy = side.energy_branch[j];
        hists.at(name + to_string(index))->Fill(energy);
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