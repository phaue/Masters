//
// Adopted from Erik - https://gitlab.au.dk/ausa/erik/e21010
// Added slight modifications

#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1I.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <algorithm>

using namespace EUtil;
using namespace std;

string input_dir  = getProjectRoot() + "/data/cal";
string output_dir = getProjectRoot() + "/data/cal";

const vector<string> dead_strips = {
    "U2F6", "U3B16", "U4F1", "U4B14", "U4B15", "U4B16"
};

int main(int argc, char* argv[]) {
  vector<string> hist_names;
  for (const string& num : {"1", "2", "3", "4", "5", "6"}) {
    for (const string& side : {"F", "B"}) {
      for (int i = 0; i < 16; i++) {
        string name = "U" + num + side + to_string(i+1);

        // skip dead strips
        if (std::find(dead_strips.begin(),dead_strips.end(),name) != dead_strips.end()) continue;

        hist_names.emplace_back(name);
      }
    }
  }
  const vector<string> p_detectors = {"P1", "P2", "P3", "P4", "P6"};
  for (const auto& p_det : p_detectors) {
      hist_names.emplace_back(p_det);
  }
  std::unique_ptr<TFile> in = make_unique<TFile>((input_dir + "/3ah.root").c_str(), "read");
  std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/3aha.root").c_str(), "recreate");

  for (const auto& hist_name : hist_names) {
    in->cd();
    TH1D *h1; in->GetObject(hist_name.c_str(), h1);

    out->cd();
    h1->Write();

    const Int_t N = h1->GetNbinsX() - 2; // first index is 'underflow', last index is 'overflow'
    Double_t spec_in[N];
    for (Int_t i = 0; i < N; i++) {
      spec_in[i] = h1->GetBinContent(i+1);
    }

    auto *s = new TSpectrum();
    Double_t spec_out[N];
    Int_t nfound = s->SearchHighRes(spec_in, spec_out, N,
                                    8, 2, kTRUE,
                                    3, kTRUE, 3);

    auto *h2 = new TH1D((hist_name + "A").c_str(), (hist_name + "A").c_str(),
                        N, h1->GetBin(0), h1->GetBin(N - 1));
    for (int i = 0; i < N; i++) h2->SetBinContent(i, spec_out[i]);
    h2->Write();

    Double_t *xpeaks = s->GetPositionX();
    Double_t pos[nfound]; Double_t amp[nfound];
    Int_t bin;
    for (int i = 0; i < nfound; i++) {
      bin = (Int_t) lround(xpeaks[i] + 0.5);
      pos[i] = h1->GetBinCenter(bin);
      amp[i] = h2->GetBinContent(bin);
    }
    auto *g1 = new TGraph(nfound, xpeaks, amp);
    g1->SetName((hist_name + "P").c_str());
    g1->SetTitle((hist_name + "P").c_str());
    g1->Write();
  }

  out->Close();
  in->Close();

  return 0;
}