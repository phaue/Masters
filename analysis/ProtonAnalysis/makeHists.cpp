#include <filesystem>
#include <vector>
#include <string>
#include <iostream>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TMarker.h>
#include "projectutil.h"
#include <TGraph.h>
#include <TArrow.h>



namespace fs = std::filesystem;
using namespace std;
using namespace EUtil;

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

int main() {
    string input_dir = getProjectRoot() + "/data/singleprotons/Al";
    string output_dir = getProjectRoot() + "/data/protonanalysis";

    // Step 1: Load matching ROOT files into TChain
    vector<string> files = getMatchingFiles(input_dir, "mlio.root");
    auto* c = new TChain("a");
    for (const auto& file : files) {
        cout << "Chaining " << file << endl;
        c->Add(file.c_str());
    }

    unsigned short id[100];
    double E[100];
    unsigned int mul = 0;

    c->SetBranchAddress("mul", &mul);
    c->SetBranchAddress("id", id);
    c->SetBranchAddress("E", E);

    // Step 3: Prepare histogram
    int nBins = 850;
    double xLow = 0.0, xHigh = 8500.0;

    auto* energyHist = new TH1D("energyHist", "Energy for id 1,2,3", nBins, xLow, xHigh);

    // Step 4: Loop over events and fill histogram
    Long64_t nEntries = c->GetEntries();
    for (Long64_t i = 0; i < c->GetEntries(); ++i) {
    c->GetEntry(i);
    for (unsigned int j = 0; j < mul; ++j) {
        if (id[j] == 1 || id[j] == 2 || id[j] == 3) {
            energyHist->Fill(E[j]);
            }
        }
    }
// --- Subhistogram Peak Search with Overlap ---
    const int nSub = 17;
    const double overlapFraction = 0.25; // 20% overlap
    const double resolution = 2.0;
    const int overlapBins = static_cast<int>(overlapFraction * (nBins / nSub));
    const int subWidthBins = nBins / nSub;
    const double binWidth = (xHigh - xLow) / nBins;
    const double peakMergeTolerance = binWidth; // merge peaks within 1 bin width

    vector<double> allPeaks;
    TSpectrum spectrum(100);
    spectrum.SetResolution(resolution);

    for (int i = 0; i < nSub; ++i) {
        int binStart = i * subWidthBins + 1;
        int binEnd = binStart + subWidthBins - 1;

        if (i > 0) binStart -= overlapBins;
        if (i < nSub - 1) binEnd += overlapBins;
        if (binStart < 1) binStart = 1;
        if (binEnd > nBins) binEnd = nBins;

        TH1D* subHist = (TH1D*)energyHist->Clone(Form("subHist_%d", i));
        subHist->SetTitle(Form("Subhistogram %d", i));

        for (int b = 1; b < binStart; ++b) subHist->SetBinContent(b, 0);
        for (int b = binEnd + 1; b <= nBins; ++b) subHist->SetBinContent(b, 0);

        int nFound = spectrum.Search(subHist);
        double* positions = spectrum.GetPositionX();

        for (int p = 0; p < nFound; ++p) {
            double pos = positions[p];
            if (pos >= subHist->GetBinLowEdge(binStart) &&
                pos <= subHist->GetBinLowEdge(binEnd) + subHist->GetBinWidth(binEnd)) {
                allPeaks.push_back(pos);
            }
        }

        delete subHist;
    }

    // --- Deduplicate peaks ---
    sort(allPeaks.begin(), allPeaks.end());
    vector<double> uniquePeaks;
    for (size_t i = 0; i < allPeaks.size(); ++i) {
        if (i == 0 || fabs(allPeaks[i] - allPeaks[i - 1]) > peakMergeTolerance) {
            uniquePeaks.push_back(allPeaks[i]);
        }
    }

    cout << "Unique peaks found across all subhistograms: " << uniquePeaks.size() << endl;
    for (size_t i = 0; i < uniquePeaks.size(); ++i) {
        cout << "Peak " << i + 1 << ": Position = " << uniquePeaks[i] << " keV" << endl;
    }

    TCanvas* c1 = new TCanvas("c1", "Energy with Peak Arrows", 800, 600);
    energyHist->Draw(); // Draw the histogram

    for (double x : uniquePeaks) {
        double y = energyHist->GetBinContent(energyHist->FindBin(x));
        auto* arrow = new TArrow(x, y + 0.02 * energyHist->GetMaximum(), x, y + 0.1 * energyHist->GetMaximum(), 0.02, "|>");
        arrow->SetLineColor(kRed);
        arrow->SetFillColor(kRed);
        arrow->SetLineWidth(2);
        arrow->Draw();
    }

    // --- Step 5: Save histograms to file ---
    TFile* outFile = new TFile((output_dir + "/PHists.root").c_str(), "RECREATE");
    energyHist->Write("energyHist");
    c1->Write("energyHistA");

    // If you want to write the TSpectrum object too:
    //spectrum.Write("spectrum");

    outFile->Close();
    std::cout << "Histograms saved to PHists.root" << std::endl;
    return 0;
}


