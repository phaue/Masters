#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1D.h>
#include <filesystem>
#include <vector>
#include <string>
#include <memory>
#include <map>
#include <TSpectrum.h>
#include <TF1.h> 
#include <fstream>
#include <numeric>


namespace fs = std::filesystem;
using namespace EUtil;
using namespace std;

vector<string> getMatchingFiles(const string& dir, const string& suffix) {
    vector<string> matched;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.is_regular_file()) {
            string filename = entry.path().filename().string();
            if (filename.size() >= suffix.size() &&
                filename.compare(filename.size() - suffix.size(), suffix.size(), suffix) == 0) {
                matched.push_back(entry.path().string());  // full path
            }
        }
    }
    return matched;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_subdir>\n";
        cerr << "Allowed values: padmatcher, padcalc\n";
        return 1;
    }

    string subdir = argv[1];
    if (subdir != "padmatcher" && subdir != "padcalc") {
        cerr << "Error: input_subdir must be either 'padmatcher' or 'padcalc'\n";
        return 1;
    }
    string input_dir = getProjectRoot() + "data/" + subdir + "/Si/";
    string output_dir = getProjectRoot() + "data/pad_calibration/";
    string output_file = output_dir + subdir + ".root";

    vector<string> hist_names;

    if (subdir == "padmatcher") {
        hist_names = {"P1ch", "P2ch", "P3ch", "P6ch"};
    } else if (subdir == "padcalc") {
        hist_names = {"P1cal", "P2cal", "P3cal", "P6cal"};
    } else {
        cerr << "Unexpected input name: " << subdir << endl;
        return 1;  // or handle error
    }
    map<string, unique_ptr<TH1D>> summed_hists;

    vector<string> files = getMatchingFiles(input_dir, ".root");
    for (const auto& file : files) {
        cout << "Reading from " << file << endl;
        unique_ptr<TFile> in = make_unique<TFile>(file.c_str(), "read");
        if (in->IsZombie()) {
            cerr << "Error opening file: " << file << endl;
            continue;
        }

        for (const auto& name : hist_names) {
            TH1D* h = nullptr;
            in->GetObject(name.c_str(), h);

            if (!h) {
                cerr << name << " not found in " << file << endl;
                continue;
            }

            if (!summed_hists[name]) {
                summed_hists[name].reset(dynamic_cast<TH1D*>(h->Clone(name.c_str())));
                summed_hists[name]->SetDirectory(nullptr);
            } else {
                summed_hists[name]->Add(h);
            }
        }
    }

    cout << "Writing summed histograms + peak search to: " << output_file << endl;
    unique_ptr<TFile> Hout = make_unique<TFile>(output_file.c_str(), "RECREATE");
    
    string peak_fits_filename = output_file.substr(0, output_file.find_last_of('.')) + "_peaks.txt";
    ofstream Pout(peak_fits_filename, ios::trunc);


    
    for (const auto& name : hist_names) {
        if (!summed_hists[name]) continue;

        TH1D* h1 = summed_hists[name].get();
        h1->Write();  // original histogram

        const Int_t N = h1->GetNbinsX() - 2;
        vector<Double_t> spec_in(N);
        for (Int_t i = 0; i < N; i++) {
            spec_in[i] = h1->GetBinContent(i + 1);
        }

        auto* s = new TSpectrum();
        vector<Double_t> spec_out(N, 0.0);
        Int_t nfound = s->SearchHighRes(spec_in.data(), spec_out.data(), N,
                                        8, 2, kTRUE, 3, kTRUE, 3);

        const int max_peaks = 4; // Set your desired maximum number of peaks
        if (nfound > max_peaks) {
            nfound = max_peaks;
        }

        // Construct filtered histogram
        auto* h2 = new TH1D((name + "A").c_str(), (name + "A").c_str(),
                            h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());
        for (Int_t i = 0; i < N; i++) {
            h2->SetBinContent(i + 1, spec_out[i]);
        }
        h2->Write();

        Double_t* xpeaks = s->GetPositionX();
        cout << "Found " << nfound << " peaks in " << name << ": ";
        for (int i = 0; i < nfound; ++i) {
            cout << xpeaks[i]*h1->GetBinWidth(1) << " ";
        }
        cout << endl;

    auto* fitHist = (TH1D*)h1->Clone((name + "_fits").c_str());
    fitHist->SetTitle((name + " with Gaussian Fits").c_str());
    fitHist->SetLineColor(kBlack);

    vector<pair<double, double>> peaks;

vector<double> xpeaks_vec(nfound);
vector<double> sigmas_vec(nfound);

for (int i = 0; i < nfound; ++i) {
    xpeaks_vec[i] = xpeaks[i]*h1->GetBinWidth(1);

    double peakX = xpeaks[i]*h1->GetBinWidth(1);
    int bin = h1->FindBin(peakX);
    double height = h1->GetBinContent(bin);
    double sigmaGuess = 15 * h1->GetBinWidth(1);

    TF1* fit = new TF1(Form("gaus_%s_%d", name.c_str(), i), "gaus", peakX - 3 * sigmaGuess, peakX + 3 * sigmaGuess);
    fit->SetParameters(height, peakX, sigmaGuess);
    fit->SetLineColor(kRed + i % 6);
    fit->SetLineWidth(2);

    fitHist->Fit(fit, "RQ+");

    sigmas_vec[i] = fit->GetParameter(2); // store fitted sigma

    delete fit;
}

// --- Step 2: Sort using indices ---
vector<size_t> indices(nfound);
iota(indices.begin(), indices.end(), 0);

sort(indices.begin(), indices.end(), [&](size_t i1, size_t i2) {
    return xpeaks_vec[i1] < xpeaks_vec[i2];
});

// --- Step 3: Output sorted results ---
for (size_t rank = 0; rank < indices.size(); ++rank) {
    size_t i = indices[rank];
    Pout << name << " " << rank + 1 << " "
         << xpeaks_vec[i] << " "
         << sigmas_vec[i] << "\n";
}

    fitHist->Write();
    delete s;
}
Pout.close();
Hout->Close();
return 0;
}