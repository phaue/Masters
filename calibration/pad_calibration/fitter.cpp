#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <ausa/calibration/Source.h>
#include <ausa/calibration/peakfinder/SimplePeakFinder.h>
#include <ausa/calibration/peakfinder/GaussianMultiPeakFinder.h>
#include "ausa/util/HistogramUtil.h"
#include <ausa/json/IO.h>
#include <ausa/util/Resource.h>





using namespace std;
using namespace AUSA::Calibration;


int main(int argc, char *argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <root_file>" << endl;
        return 1;
    }

    // Open the ROOT file containing the histograms
    TFile *file = TFile::Open(argv[1]);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << argv[1] << endl;
        return 1;
    }


TH1D *histP1 = (TH1D*)file->Get("P1ch");
TH1D *histP2 = (TH1D*)file->Get("P2ch");
TH1D *histP3 = (TH1D*)file->Get("P3ch");

if (!histP1 || histP1->GetEntries() == 0 || 
    !histP2 || histP2->GetEntries() == 0 || 
    !histP3 || histP3->GetEntries() == 0) {
    
    // Extract alternative histograms from the file
    histP1 = (TH1D*)file->Get("P1cal");
    histP2 = (TH1D*)file->Get("P2cal");
    histP3 = (TH1D*)file->Get("P3cal");

    if (!histP1 || histP1->GetEntries() == 0 || 
        !histP2 || histP2->GetEntries() == 0 || 
        !histP3 || histP3->GetEntries() == 0) {
        std::cerr << "Error extracting histograms!" << std::endl;
        return 1;
    }
}



    Source source{  
        { {2100}, {3100},  {6000}/*, {6000}, {7500}, {7600}*/}, // 4 peaks 
        TVector3(0,0,0),                       // Position
        TVector3(0,0,-1),                       // Direction
        1,                                   // Thickness
        1                                  // Transverse size
      };
      //auto source = AUSA::JSON::readSourceFromJSON("/home/haue/repositories/Masters/analysis/padcal/source.json");
/*
      //auto simple = std::make_unique<SimplePeakFinder>(200);
      SimplePeakFinder simple(200);
      //auto finder = std::make_unique<GaussianMultiPeakFinder>(move(simple));
      auto peaksP1 = simple.findPeaks(*histP1, source);
      auto peaksP2 = simple.findPeaks(*histP2, source);
      auto peaksP3 = simple.findPeaks(*histP3, source);*/

        auto simple = std::make_unique<SimplePeakFinder>(600);
        auto finder = std::make_unique<GaussianMultiPeakFinder>(move(simple));
        finder -> setRangeMultiplier(0.8);
        auto peaksP1 = finder->findPeaks(*histP1, source);
        auto peaksP2 = finder->findPeaks(*histP2, source);
        auto peaksP3 = finder->findPeaks(*histP3, source);
      

      

        cout << "# Peaks found from fitting routine for the three pads P1,P2 and P3" << endl;
      for (size_t i =0; i<peaksP1.peakCount(); i++){
        cout << "\t" << peaksP1.getPosition(i);
     }
     cout  << "\t" << "#Mean P1"<< endl;

     for (size_t i =0; i<peaksP1.peakCount(); i++){
        cout << "\t" << peaksP1.getUncertaincy(i);
     }
     cout  << "\t" << "#Mean Error P1"<< endl;

     for (size_t i =0; i<peaksP1.peakCount(); i++){
      cout << "\t" << AUSA::estimateSigma(*histP1, peaksP1.getPosition(i));
      }
      cout  << "\t" << "#Sigma P1"<< endl;


      for (size_t i =0; i<peaksP1.peakCount(); i++){
        cout << "\t" << peaksP2.getPosition(i);
     }
     cout  << "\t" << "#Mean P2"<< endl;

      for (size_t i =0; i<peaksP1.peakCount(); i++){
        cout << "\t" << peaksP2.getUncertaincy(i);
     }
     cout  << "\t" << "#Mean Error P2"<< endl;

     for (size_t i =0; i<peaksP1.peakCount(); i++){
      cout << "\t" << AUSA::estimateSigma(*histP2, peaksP2.getPosition(i));
      }
      cout  << "\t" << "#Sigma P2"<< endl;

      for (size_t i =0; i<peaksP3.peakCount(); i++){
        cout << "\t" << peaksP3.getPosition(i);
     }
     cout  << "\t" << "#Mean P3"<< endl;

      for (size_t i =0; i<peaksP3.peakCount(); i++){
        cout << "\t" << peaksP3.getUncertaincy(i);
     }
     cout  << "\t" << "#Mean Error P1"<< endl;

     for (size_t i =0; i<peaksP3.peakCount(); i++){
      cout << "\t" << AUSA::estimateSigma(*histP3, peaksP3.getPosition(i));
      }
      cout  << "\t" << "#Sigma P3"<< endl;


    TFile *outputFile = TFile::Open("fit_results.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        cerr << "Error creating output file!" << endl;
        return 1;
    }
    
    outputFile->cd();
    histP1->Write();
    histP2->Write();
    histP3->Write();
    outputFile->Close();
    
    file->Close();

    return 0;
}
