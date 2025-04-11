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
          else continue;
        }//if dsssd and pad is partnered
      }//for each pad hit in telescope backside candidates
    }// for each dsssd hit in telescope frontside candidates
  }//specificanalysis
};//CalibrationMatcher

class BackwardCalculation : public CalibrationAnalysis {
  public:
  BackwardCalculation(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
                    string _isotopetype, bool _events_matched)
        : CalibrationAnalysis(_setupSpecs, _target, output, _isotopetype,_events_matched), iso(_isotopetype){}
  
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
            for (size_t i =0;i<peaks1_means.size(); i++){
              if(pad_hit->Edep > peaks1_means[i]-2*peaks1_sig[i] && pad_hit->Edep < peaks1_means[i]+2*peaks1_sig[i]){
                bool success = treatBackwardHit(dsssd_hit, pad_hit ,real_peaks[i]);
                if (success){
                  histP1->Fill(pad_hit->Ecal);
                  addTelescopeHit(dsssd_hit,pad_hit);
                  
                }

              }
              else continue;
            }
          }//if P2
          if(pad_det->getName()=="P2"){
            for (size_t i =0;i<peaks2_means.size(); i++){
              if(pad_hit->Edep > peaks2_means[i]-2*peaks2_sig[i] && pad_hit->Edep < peaks2_means[i]+2*peaks2_sig[i]){
                bool success = treatBackwardHit(dsssd_hit, pad_hit ,real_peaks[i]);
                if (success){
                  histP2->Fill(pad_hit->Ecal);
                  addTelescopeHit(dsssd_hit,pad_hit);
                  
                }

              }
              else continue;
            }
          }//if P3
          if(pad_det->getName()=="P3"){
            for (size_t i =0;i<peaks3_means.size(); i++){
              if(pad_hit->Edep > peaks3_means[i]-2*peaks3_sig[i] && pad_hit->Edep < peaks3_means[i]+2*peaks3_sig[i]){
                bool success = treatBackwardHit(dsssd_hit, pad_hit ,real_peaks[i]);
                if (success){
                  histP3->Fill(pad_hit->Ecal);
                  addTelescopeHit(dsssd_hit,pad_hit);
                  
                }

              }
              else continue;
            }
          }//if P3

        }//if dsssd and pad is partnered
      }//for each pad hit in telescope backside candidates
    }// for each dsssd hit in telescope frontside candidates
  
        ///These peaks are for Mg
        //They are found by running the CalibrationMatcher on pad-dummy calibrated data and running the fitter
        //A dynamic read of the file would be a great addition such that it is easy to make changes
if (iso=="Mg"){
  peaks1_means= {895.567, 1175.43, 1677.54};
  peaks1_sig ={41.0669 ,26.9267, 20.0444};
  peaks2_means= {879.403, 1223.8,  1815.38};
  peaks2_sig ={49.9033 ,27.3036 ,12.8827};
  peaks3_means= {943.037 ,1259.37, 1821.78};
  peaks3_sig ={45.537 , 25.8725 ,12.2788};
  real_peaks ={3843, 4675, 6231};
  }
else if (iso=="Si"){
  peaks1_means ={702.102, 981.697, 1177.82, 1422.2};
  peaks1_sig ={41.6016, 26.1971, 28.3199, 20.7001};
  peaks2_means ={641.157, 987.765, 1223.12, 1511.87};
  peaks2_sig ={49.657,  31.2646, 26.6186, 24.3743};
  peaks3_means = {724.738, 1039.4,  1257.86, 1531.15};
  peaks3_sig = {42.2379, 22.9041, 29.3639, 21.6494};
  real_peaks = {3338, 4089, 4652, 5403};
  }
else {
  cout << "How did you get this far? the program should terminate if no isotopetype is specifed..." << endl;
  }

}//specificanalysis

string iso;
vector<double> peaks1_means, peaks1_sig,peaks2_means, peaks2_sig, peaks3_means,peaks3_sig, real_peaks;

};//BackwardsCalculation



////should read in from cfg file, there it should be specified what isotope it is and the path to directories
//and not least which analysis we want to carry out.

int main(int argc, char *argv[]) {


  cout << "Config file path: " << getProjectRoot() + "analysis/padcal/" + getBasename(argv[1]) << endl;
  prepareFileIO(getProjectRoot() + "analysis/padcal/" + getBasename(argv[1]));
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
    analysis = make_shared<BackwardCalculation>(setup, target, &output, isotopetype, events_matched);
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

