
#include <string>
#include <vector>
#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1D.h>
#include <TSpectrum.h>
#include <TGraph.h>
#include <algorithm>
#include <TF1.h>
#include <filesystem> 
#include <memory>

using namespace EUtil;
using namespace std;
using namespace std::filesystem;

vector<string> find_files_with_ending(const string& directoryPath, const string& fileEnding) {
    vector<string> matching_files;

    // Check if the provided path actually exists and is a directory
    if (!exists(directoryPath) || !is_directory(directoryPath)) {
        cerr << "Error: Directory does not exist or is not a directory: " << directoryPath << endl;
        return matching_files; // Return empty vector
    }

    // Use a directory_iterator to loop through all entries in the folder
    for (const auto& entry : directory_iterator(directoryPath)) {
        // We are only interested in regular files, not subdirectories
        if (entry.is_regular_file()) {
            string filename = entry.path().filename().string();

            // Check if the filename ends with the desired pattern
            // A simple way to check for "*mlio.root" is to see if it ends with "mlio.root"
            if (filename.length() >= fileEnding.length() && 
                filename.compare(filename.length() - fileEnding.length(), fileEnding.length(), fileEnding) == 0) 
            {
                // If it matches, add the full path to our vector
                matching_files.push_back(entry.path().string());
            }
        }
    }

    return matching_files;
}
string input_dir  = getProjectRoot() + "/data/singleprotons/Al";
string output_dir = getProjectRoot() + "/data/singleprotons/";

vector<string> files = find_files_with_ending(input_dir, "mlio.root");

int main(int argc, char* argv[]) {

    std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/peakchecks.root").c_str(), "recreate");


    auto *c = new TChain("a");
    for (const auto& file : files) {
        cout << "Chaining " << file << endl;
    c->Add((file).c_str());
  }


    TH1D* hist = new TH1D("E", "E", 700, 1000, 8000);


    const int MAX_MULTIPLICITY = 100; 
    UInt_t mul; UShort_t id[MAX_MULTIPLICITY];
    double_t E[MAX_MULTIPLICITY];

    // Set the branch addresses for the TChain
    c->SetBranchAddress("mul", &mul);
    c->SetBranchAddress("E", E);
    c->SetBranchAddress("id", id);

    // Loop over all events in the chain
    Long64_t n_entries = c->GetEntries();
    cout << "Processing " << n_entries << " total events..." << endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        c->GetEntry(i);
        // Loop over particles within the current event
        for (int j = 0; j < mul; ++j) {
            // Apply your cut: if id is 3, fill the histogram with E
            if (id[j] == 3 ) {
                hist->Fill(E[j]);
            }
        }
    }
    
    cout << "Finished processing events." << endl;  


    out->cd();
    hist->Write();

    const Int_t N = hist->GetNbinsX(); // first index is 'underflow', last index is 'overflow'
    Double_t spec_in[N];
    for (Int_t i = 0; i < N; i++) {
      spec_in[i] = hist->GetBinContent(i+1);
    }

    auto *s = new TSpectrum();
    Double_t spec_out[N];
    Int_t nfound = s->SearchHighRes(spec_in, spec_out, N,
                                    1, 1, kTRUE,
                                    3, kTRUE, 1);

    auto *hist2 = new TH1D("A", "A",
                        N, 1000, 8000);
    for (int i = 0; i < N; i++) hist2->SetBinContent(i, spec_out[i]);
    hist2->Write();


    Double_t *xpeaks = s->GetPositionX();
    Double_t amp[nfound];
    Int_t bin;
    Double_t real_x_positions[nfound];

    for (int i = 0; i < nfound; i++) {
        
        double peak_bin_number = xpeaks[i];

        bin = static_cast<Int_t>(peak_bin_number + 1.0);

        real_x_positions[i] = hist->GetXaxis()->GetBinCenter(bin);

        amp[i] = hist->GetBinContent(bin); 
    }

auto *g1 = new TGraph(nfound, real_x_positions, amp);

g1->SetName("P");
g1->SetTitle("P");
g1->Write();


  out->Close();
  return 0;
}