//
// Created by erik on 11/12/23.
//

#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1I.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <algorithm>
#include <TF1.h>

using namespace EUtil;
using namespace std;

string input_dir  = getProjectRoot() + "/data/cal";
string output_dir = getProjectRoot() + "/data/cal";

const vector<string> dead_strips = {
    "U5F2"
};

int main(int argc, char* argv[]) {
  vector<string> hist_names;
  for (const string& num : {"4"}) {
    for (const string& side : {"F", "B"}) {
      for (int i = 0; i < 16; i++) {
        string name = "U" + num + side + to_string(i+1);

        // skip dead strips
        if (std::find(dead_strips.begin(),dead_strips.end(),name) != dead_strips.end()) continue;

        hist_names.emplace_back(name);
      }
    }
  }

  std::unique_ptr<TFile> in = make_unique<TFile>((input_dir + "/ph.root").c_str(), "read");
  std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/pha.root").c_str(), "recreate");

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
    Double_t pos[nfound]; Double_t amp[nfound]; double_t poserr[nfound]; double_t fitpos[nfound];
    Int_t bin;
    for (int i = 0; i < nfound; i++) {
      bin = (Int_t) lround(xpeaks[i] + 0.5);
      pos[i] = h1->GetBinCenter(bin); // gives the approximate position of the peaks
      amp[i] = h2->GetBinContent(bin); //gives the relative amplitude of the peaks -->> thats why we use the "cleaned" spectrum since we filter out noise to get a better sense of the relative sizes

      // Perform Gaussian fit around the peak
      double range = h1->GetBinWidth(1) * 5; // Fit range (e.g., 3 bins around the peak)
      TF1 *gausFit = new TF1("gausFit", "gaus", h1->GetBinCenter(bin) - range, h1->GetBinCenter(bin) + range);
      h1->Fit(gausFit, "RQ"); // Quiet fit
      poserr[i] = gausFit->GetParError(1); // Uncertainty on the mean
      fitpos[i] = gausFit->GetParameter(1);

    }
    auto *g1 = new TGraph(nfound, xpeaks, amp);
    auto *g2 = new TGraph(nfound, fitpos, poserr);
    g1->SetName((hist_name + "P").c_str());
    g1->SetTitle((hist_name + "P").c_str());
    g2->SetName((hist_name + "F").c_str());
    g2->SetTitle((hist_name + "F").c_str());
    g1->Write();
    g2->Write();
  }

  out->Close();
  in->Close();

  return 0;
}