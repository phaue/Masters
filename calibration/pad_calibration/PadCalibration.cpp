#include <libconfig.h++>
#include "projectutil.h"
#include "Hit.h"
#include "AnalysisConfig.h"
#include "CalibrationAnalysis.h"
#include <gsl/gsl_vector.h>
#include "fstream"

#include <TROOT.h>
#include <Math/Interpolator.h>
#include <ausa/json/IO.h>
#include <ausa/setup/Target.h>
#include <ausa/eloss/Default.h>
#include <ausa/util/FileUtil.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/calibration/Source.h>
#include <ausa/calibration/peakfinder/SimplePeakFinder.h>
#include <ausa/calibration/peakfinder/GaussianMultiPeakFinder.h>

using namespace std;
using namespace libconfig;
using namespace EUtil;
using namespace Detectors_;
using namespace AUSA;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace ROOT::Math;
using namespace ROOT::Math::Interpolation;
using namespace AUSA::Calibration;


class CalibrationMatcher : public CalibrationAnalysis {
  public:
  CalibrationMatcher(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
                    string _isotopetype, bool _events_matched)
        : CalibrationAnalysis(_setupSpecs, _target, output, _isotopetype,_events_matched){}

  void specificAnalysis() override {
    if (hits.empty()) return;

    unordered_set<Hit*> telescope_frontside_candidates;
    unordered_set<Hit*> telescope_backside_candidates;

    for(auto &hit : hits) {
      auto det = hit.detector;
      switch(det->getType()) {
        case DSSSD:
          if(det->hasPartner()){
            telescope_frontside_candidates.emplace(&hit);
            }
            else {
              continue;
            }
              break;
        case Pad:
          telescope_backside_candidates.emplace(&hit);
          break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
      }//switch
    }//for hit in hits

    for(auto dsssd_hit : telescope_frontside_candidates) {
      auto dsssd_det = dsssd_hit->detector;
      for(auto pad_hit : telescope_backside_candidates) {
        auto pad_det = pad_hit->detector;
        if(dsssd_det->getPartner() == pad_det && pad_hit->Edep>500){
          if(pad_det->getName()=="P1"){
            //cout << "P1 with Ech>500" << endl;
              histP1->Fill(pad_hit->Edep);
              addPadHit(pad_hit);  
            }
          else if (pad_det->getName()=="P2"){
             histP2->Fill(pad_hit->Edep);
             addPadHit(pad_hit);}
          else if (pad_det->getName()=="P3"){
             histP3->Fill(pad_hit->Edep);
             addPadHit(pad_hit);}
          else if (pad_det->getName()=="P6"){
             histP6->Fill(pad_hit->Edep);
             addPadHit(pad_hit);}
          else continue;
        }//if dsssd and pad is partnered
      }//for each pad hit in telescope backside candidates
    }// for each dsssd hit in telescope frontside candidates
  }//specificanalysis
};//CalibrationMatcher

class BackwardCalculation : public CalibrationAnalysis {
  public:
  BackwardCalculation(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
                    string _isotopetype, bool _events_matched, const string& _peak_path)
        : CalibrationAnalysis(_setupSpecs, _target, output, _isotopetype,_events_matched), peak_path(_peak_path)
        {
    ifstream infile(peak_path);
    if (!infile.is_open()) {
        cerr << "Could not open peak file: " << peak_path << endl;
        return;
    }

    string detector;
    int ch;
    double peak, sigma, fitmean;

    while (infile >> detector >> ch >> peak >> fitmean >> sigma) {
        string detID = detector.substr(0, 2); // Get "P1", "P2", etc.
        detectorPeaks[detID].means.push_back(fitmean);
        detectorPeaks[detID].sigmas.push_back(sigma);
        //cout << "detector " << detector << "   ch " << ch << "    peak " << fitmean << "    sigma " << sigma << endl;  
    }
    }

  
    void specificAnalysis() override {
    if (hits.empty()) return;

    unordered_set<Hit*> telescope_frontside_candidates;
    unordered_set<Hit*> telescope_backside_candidates;

    for(auto &hit : hits) {
      auto det = hit.detector;
      switch(det->getType()) {
        case DSSSD:
          if(det->hasPartner()){
            telescope_frontside_candidates.emplace(&hit);
            }
            else {
              continue;
            }
              break;
        case Pad:
          telescope_backside_candidates.emplace(&hit);
          break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
      }//switch
    }//for hit in hits

for (auto dsssd_hit : telescope_frontside_candidates) {
    auto dsssd_det = dsssd_hit->detector;

    for (auto pad_hit : telescope_backside_candidates) {
        auto pad_det = pad_hit->detector;

        // 1) Check partner and energy threshold
        if (dsssd_det->getPartner() != pad_det || pad_hit->Edep <= 500){
            continue;}

        // 2) Look up “P1”, “P2” or “P3” in your map
        string detName = pad_det->getName();
        auto it = detectorPeaks.find(detName);
        if (it == detectorPeaks.end())
            continue;

        // 3) Reference that detector’s peak data
        PeakData &peaks = it->second;

        // 4) Loop over all peaks for this detector
        for (size_t i = 0; i < peaks.means.size(); ++i) {
            double mean = peaks.means[i];
            double sig  = peaks.sigmas[i];
            //out << "mean" << mean << "    sig" << sig << "    realpeak " << real_peaks[i] << endl;
            if (pad_hit->Edep > mean - 1.*sig &&
                pad_hit->Edep < mean + 1.*sig)
            {
                bool success = treatBackwardHit(dsssd_hit, pad_hit, real_peaks[i]);
                if (success) {
                    // 5) Fill the correct histogram based on detName
                    if (detName == "P1")      histP1->Fill(pad_hit->Ecal);
                    else if (detName == "P2") histP2->Fill(pad_hit->Ecal);
                    else if (detName == "P3") histP3->Fill(pad_hit->Ecal);
                    else if (detName == "P6") histP6->Fill(pad_hit->Ecal);

                    addTelescopeHit(dsssd_hit, pad_hit);
                }
            }
        }
    }
}
        }
private:
  struct PeakData {
    vector<double> means;
    vector<double> sigmas;
  };

  map<string, PeakData> detectorPeaks;
  string peak_path;
  const vector<double> real_peaks = {
    3337.75,
    4089.18,
    4651.19,
    5402.61
  };
};



////should read in from cfg file, there it should be specified what isotope it is and the path to directories
//and not least which analysis we want to carry out.

int main(int argc, char *argv[]) {


  cout << "Config file path: " << getProjectRoot() + "calibration/pad_calibration/" + getBasename(argv[1]) << endl;
  prepareFileIO(getProjectRoot() + "calibration/pad_calibration/" + getBasename(argv[1]));
  auto setup = JSON::readSetupFromJSON(setup_path);
  auto target = make_shared<Target>(JSON::readTargetFromJSON(target_path));

  system(("mkdir -p " + output_path_dir).c_str());

  vector<string> input;
  //auto setup = JSON::readSetupFromJSON("/home/haue/repositories/Masters/setup/cal_setup.json");
  //auto target = make_shared<Target>(JSON::readTargetFromJSON("/home/haue/repositories/Masters/setup/target.json"));

  for (int i = 2; i < argc; i++) {
    //check if argv is a file
    ifstream file(argv[i]);
    if (!file.good()) {
        cerr << "File " << argv[i] << " not found!" << endl;
        return 1;
    }
    cout << "using file " << argv[i] << endl;
    input.push_back(argv[i]);
  }

for (auto &in : input) {
  clock_t start = clock();
  SortedReader reader{*setup};
  reader.add(in);
  reader.setVerbose(true);
  auto base = stripFileExtension(extractFileName(in));
  string stem = getStem(in);
  TString outfile = (output_path_dir + "/" + stem + "a.root").c_str(); //specifies the name of the output
  TFile output(outfile, "RECREATE");
  prepareAnalysis();
  shared_ptr<CalibrationAnalysis> analysis;
  if (specificAnalysis == "CalibrationMatcher"){
    analysis = make_shared<CalibrationMatcher>(setup, target, &output, isotopetype, events_matched);
  }
  else if (specificAnalysis == "BackwardCalculation"){
    if (peak_path.empty()){
      cerr << "No peak path run CalibrationMatcher first and then finder before continuing" << endl;
    }
    analysis = make_shared<BackwardCalculation>(setup, target, &output, isotopetype, events_matched, peak_path);
  }
  else {
    cerr << "Type of analysis not recognizable -- recheck config file -- Aborting analysis." << endl;
    abort();
  }

  printConfig();
  cout << endl << endl << "Reading from: " << in << ' ' << endl;
  cout << "Printing to:  " << outfile << endl;
  reader.attach(analysis);
  reader.run();

  //print the time it took to analyze the file.
  clock_t stop = clock();
  double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
  printf("\nTime elapsed: %.5f\n", elapsed);
}

return 0; // SUCCESS!
}

