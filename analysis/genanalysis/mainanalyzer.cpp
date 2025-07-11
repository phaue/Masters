#include <libconfig.h++>
#include "projectutil.h"
#include "Hit.h"
#include "AnalysisConfig.h"
#include "GeneralAnalysis.h"
#include "U1analysis.h"
#include <gsl/gsl_vector.h>
#include <TROOT.h>
#include <Math/Interpolator.h>
#include <ausa/json/IO.h>
#include <ausa/setup/Target.h>
#include <ausa/eloss/Default.h>
#include <ausa/util/FileUtil.h>
#include <ausa/sort/SortedReader.h>
using namespace std;
using namespace libconfig;
using namespace EUtil;
using namespace Detectors_;
using namespace AUSA;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace ROOT::Math;
using namespace ROOT::Math::Interpolation;

class BananaMaker : public GeneralAnalysis {
  public:
  BananaMaker(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output, string _isotopetype,
                    bool _exclude_hpges = false, bool _include_DSSSD_rim = false, bool _include_spurious_zone = false,
                    bool _include_banana_cuts=false, bool _include_beta_region =false)
        : GeneralAnalysis(_setupSpecs, _target, output, _isotopetype, _exclude_hpges,
                          _include_DSSSD_rim, _include_spurious_zone, _include_banana_cuts, _include_beta_region) {}

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
            else { //U5 doesnt have a partner so we need to treat it aswell
              treatDSSSDHit(&hit);
              addDSSSDHit(&hit);
              }
              break;
        case Pad:
          telescope_backside_candidates.emplace(&hit);
          break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
      }//switch
    }//for hit in hits

    unordered_set<Hit*> telescope_frontside_successes;
    unordered_set<Hit*> telescope_backside_successes;
    for(auto dsssd_hit : telescope_frontside_candidates) {
      auto dsssd_det = dsssd_hit->detector;
      for(auto pad_hit : telescope_backside_candidates) {
        auto pad_det = pad_hit->detector;
        if(dsssd_det->getPartner() == pad_det){
          bool telescope_success = treatTelescopeHit(dsssd_hit, pad_hit);
          if(telescope_success){
            telescope_frontside_successes.emplace(dsssd_hit);
            telescope_backside_successes.emplace(pad_hit);
            addTelescopeHit(dsssd_hit, pad_hit);
/*
            Each hit in the DSSSD is paired with all pad hits and this is done for each hit in the DSSSD
            This creates multiple pairings between the two
 */
          }//if telescope hit is a success
        }//if dsssd and pad is partnered
      }//for each pad hit in telescope backside candidates
    }// for each dsssd hit in telescope frontside candidates

    for(auto hit : telescope_frontside_successes){
      telescope_frontside_candidates.erase(hit);
      }//removes all hits that were a success from the frontside candidates
  
    for(auto hit : telescope_frontside_candidates){
      treatDSSSDHit(hit);
      addDSSSDHit(hit);
      }//treats the leftover non-matched dsssd hits.
  }//specificanalysis
};//BananaMaker


class SingleProton : public GeneralAnalysis{
  public:
  SingleProton(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output, string _isotopetype,
                    bool _exclude_hpges = false, bool _include_DSSSD_rim = false, bool _include_spurious_zone = false,
                    bool _include_banana_cuts=false, bool _include_beta_region =false)
        : GeneralAnalysis(_setupSpecs, _target, output, _isotopetype, _exclude_hpges,
                          _include_DSSSD_rim, _include_spurious_zone, _include_banana_cuts, _include_beta_region) {}

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
            else { //U5 doesnt have a partner so we need to treat it aswell
              treatDSSSDHit(&hit);
              addDSSSDHit(&hit);
              }
              break;
        case Pad:
          telescope_backside_candidates.emplace(&hit);
          break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
      }//switch
    }//for hit in hits


    unordered_set<Hit*> telescope_frontside_successes;
    unordered_set<Hit*> telescope_backside_successes;
    for(auto dsssd_hit : telescope_frontside_candidates) {
      auto dsssd_det = dsssd_hit->detector;
      for(auto pad_hit : telescope_backside_candidates) {
        auto pad_det = pad_hit->detector;
        if(dsssd_det->getPartner() == pad_det){
          bool telescope_success = treatTelescopeHit(dsssd_hit, pad_hit);
          if(telescope_success){

            if(dsssd_det->getBananaCut()->isSatisfied(pad_hit->Edep, dsssd_hit->Edep)){

            
            telescope_frontside_successes.emplace(dsssd_hit);
            telescope_backside_successes.emplace(pad_hit);
            addTelescopeHit(dsssd_hit, pad_hit);
            }
          }//if telescope hit is a success
        }//if dsssd and pad is partnered
      }//for each pad hit in telescope backside candidates
    }// for each dsssd hit in telescope frontside candidates

    for(auto hit : telescope_frontside_successes){
      telescope_frontside_candidates.erase(hit);
      }//removes all hits that were a success from the frontside candidates

    for(auto hit : telescope_frontside_candidates){
      treatDSSSDHit(hit);
      addDSSSDHit(hit);
      }//treats the leftover non-matched dsssd hits.
  }//specificanalysis
};//singleproton

class GammaSpec : public GeneralAnalysis {
  public:
  GammaSpec(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output, string _isotopetype,
                    bool _exclude_hpges = false, bool _include_DSSSD_rim = false, bool _include_spurious_zone = false,
                    bool _include_banana_cuts=false, bool _include_beta_region =false)
        : GeneralAnalysis(_setupSpecs, _target, output, _isotopetype, _exclude_hpges,
                          _include_DSSSD_rim, _include_spurious_zone, _include_banana_cuts, _include_beta_region) {}

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
            else { //U5 doesnt have a partner so we need to treat it aswell
              treatDSSSDHit(&hit);
              addDSSSDHit(&hit);
              }
              break;
        case Pad:
          telescope_backside_candidates.emplace(&hit);
          break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
      }//switch
    }//for hit in hits

    unordered_set<Hit*> telescope_frontside_successes;
    unordered_set<Hit*> telescope_backside_successes;
    for(auto dsssd_hit : telescope_frontside_candidates) {
      auto dsssd_det = dsssd_hit->detector;
      for(auto pad_hit : telescope_backside_candidates) {
        auto pad_det = pad_hit->detector;
        if(dsssd_det->getPartner() == pad_det){
          bool telescope_success = treatTelescopeHit(dsssd_hit, pad_hit);
          if(telescope_success){
            telescope_frontside_successes.emplace(dsssd_hit);
            telescope_backside_successes.emplace(pad_hit);
            addTelescopeHit(dsssd_hit, pad_hit);
/*
            Each hit in the DSSSD is paired with all pad hits and this is done for each hit in the DSSSD
            This creates multiple pairings between the two
 */
          }//if telescope hit is a success
        }//if dsssd and pad is partnered
      }//for each pad hit in telescope backside candidates
    }// for each dsssd hit in telescope frontside candidates

    for(auto hit : telescope_frontside_successes){
      telescope_frontside_candidates.erase(hit);
      }//removes all hits that were a success from the frontside candidates
      
    for(auto hit : telescope_frontside_candidates){
      treatDSSSDHit(hit);
      addDSSSDHit(hit);
      }//treats the leftover non-matched dsssd hits.
  }//specificanalysis
};//Gammaspec




class AboveBananaAnalysis : public U1analysis{
  public : AboveBananaAnalysis(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
    string _isotopetype, bool _Only_U1)
:U1analysis(_setupSpecs, _target, output, _isotopetype, _Only_U1) {}

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
        else { //U5 doesnt have a partner so we need to treat it aswell
          break;
        }
          break;
      case Pad:
        telescope_backside_candidates.emplace(&hit);
        break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
  }//switch
  }//for hit in hits
  unordered_set<Hit*> telescope_frontside_successes;
  unordered_set<Hit*> telescope_backside_successes;

  for(auto dsssd_hit : telescope_frontside_candidates) {
    auto dsssd_det = dsssd_hit->detector;
    for(auto pad_hit : telescope_backside_candidates) {
      auto pad_det = pad_hit->detector;
      if(dsssd_det->getPartner() == pad_det){
        if(dsssd_det->getBananaCut()->isSatisfied(pad_hit->Edep, dsssd_hit->Edep)){            
          if(pad_hit->Edep+dsssd_hit->Edep <4100 && pad_hit->Edep+dsssd_hit->Edep > 3900){
            bool telescope_success = specialTelescopeTreatment(dsssd_hit, pad_hit, 4100);
            if(telescope_success){
              telescope_frontside_successes.emplace(dsssd_hit);
              telescope_backside_successes.emplace(pad_hit);
              addTelescopeHit(dsssd_hit, pad_hit);
              }
          }
          else if(pad_hit->Edep+dsssd_hit->Edep <5500 && pad_hit->Edep+dsssd_hit->Edep > 5200){
            bool telescope_success = specialTelescopeTreatment(dsssd_hit, pad_hit, 5400);
            if(telescope_success){
              telescope_frontside_successes.emplace(dsssd_hit);
              telescope_backside_successes.emplace(pad_hit);
              addTelescopeHit(dsssd_hit, pad_hit);
              }
          }
          else if(pad_hit->Edep+dsssd_hit->Edep <4700 && pad_hit->Edep+dsssd_hit->Edep > 4400){
            bool telescope_success = specialTelescopeTreatment(dsssd_hit, pad_hit, 4700);
            if(telescope_success){
              telescope_frontside_successes.emplace(dsssd_hit);
              telescope_backside_successes.emplace(pad_hit);
              addTelescopeHit(dsssd_hit, pad_hit);
              }
          }
          else{
            bool telescope_success = treatTelescopeHit(dsssd_hit, pad_hit);
            if(telescope_success){
              telescope_frontside_successes.emplace(dsssd_hit);
              telescope_backside_successes.emplace(pad_hit);
              addTelescopeHit(dsssd_hit, pad_hit);
              }
          }
        }//if telescope hit is a success
      }//if dsssd and pad is partnered
    }//for each pad hit in telescope backside candidates
  }  
}//specificanalysis
};//AboveBananaAnalysis

class Bananaexplorer : public U1analysis{
  public : Bananaexplorer(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
    string _isotopetype, bool _Only_U1, vector<vector<int>> _peakSets)
:U1analysis(_setupSpecs, _target, output, _isotopetype, _Only_U1), peakSets(_peakSets) {}

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
        else { //U5 doesnt have a partner so we need to treat it aswell
          break;
        }
          break;
      case Pad:
        telescope_backside_candidates.emplace(&hit);
        break;
      default: //treats the rest of the cases such as NoType and HPGe
        break;
  }//switch
  }//for hit in hits
  unordered_set<Hit*> telescope_frontside_successes;
  unordered_set<Hit*> telescope_backside_successes;

  for(auto dsssd_hit : telescope_frontside_candidates) {
    auto dsssd_det = dsssd_hit->detector;
    for(auto pad_hit : telescope_backside_candidates) {
      auto pad_det = pad_hit->detector;
      if(dsssd_det->getPartner() == pad_det){
        for (const auto& peak: peakSets){
          int lowE = peak[0];
          int upE = peak[1];
          int iniE = peak[2];
          if(pad_hit->Edep+dsssd_hit->Edep <= upE && pad_hit->Edep+dsssd_hit->Edep >= lowE){
            bool telescope_success = specialTelescopeTreatment(dsssd_hit, pad_hit, iniE);
            if(telescope_success){
              telescope_frontside_successes.emplace(dsssd_hit);
              telescope_backside_successes.emplace(pad_hit);
              addTelescopeHit(dsssd_hit, pad_hit);
              } // if telescope success
          }//if padhit
          else{
            continue;
          }}
      }//if dsssd and pad is partnered
    }//for each pad hit in telescope backside candidates
  }//for each dsssd hit
}//specificanalysis
private:
    vector<vector<int>> peakSets;
};//Bananaexplorer

class TwoProton : public GeneralAnalysis {
public:
  TwoProton(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output, string _isotopetype,
                    bool _exclude_hpges = false, bool _include_DSSSD_rim = false, bool _include_spurious_zone = false,
                    bool _include_banana_cuts=false, bool _include_beta_region =false)
        : GeneralAnalysis(_setupSpecs, _target, output, _isotopetype, _exclude_hpges,
                          _include_DSSSD_rim, _include_spurious_zone, _include_banana_cuts, _include_beta_region) {}

  void specificAnalysis() override {
    if (hits.empty()) return;

    int proton_events = 0;
    
    //finds the telescope proton events from the hits
    pair<vector<Hit*>, vector<Hit*>> telescope_candidates;
    for (auto &hit : hits) {
      auto det = hit.detector;
      switch (det->getType()) {
        case DSSSD:
          if (det->hasPartner()
              //&& det->getId() < 4 && det->getId() > 0// look only for b2p in the best detector telescopes
          ) {
            telescope_candidates.first.emplace_back(&hit);
          }
          break;
        case Pad:
          if (det->hasPartner()
              //&& det->getId() < 10 && det->getId() > 6 // look only for b2p in the best detector telescopes
          ) {
            telescope_candidates.second.emplace_back(&hit);
          }
          break;
        default:
          // if case is 'HPGe', the scalar Eg1 and/or Eg2 already contains the gamma energy and the boolean 'g' is set to 'true'
          // do nothing
          break;
      }
    }

    pair<vector<Hit*>, vector<Hit*>> telescope_successes;
    for (auto front_hit : telescope_candidates.first) {
      if (std::find(telescope_successes.first.begin(),
                    telescope_successes.first.end(), 
                    front_hit) != telescope_successes.first.end()) continue;
      auto front_det = front_hit->detector;
      for (auto back_hit : telescope_candidates.second) {
        if (std::find(telescope_successes.second.begin(),
                      telescope_successes.second.end(),
                      back_hit) != telescope_successes.second.end()) continue;
        auto back_det = back_hit->detector;
        if (front_det->getPartner() == back_det) {
          bool success = treatTelescopeHit(front_hit, back_hit);
          if (success) {
            if (front_det->getBananaCut()->isInside(back_hit->Edep, front_hit->Edep)) {
              proton_events++;
              telescope_successes.first.emplace_back(front_hit);
              telescope_successes.second.emplace_back(back_hit);
            }
          }
        }
      }
    }
    //find the deltaE contained proton events not treated from the previous telescope treatment
    unordered_set<Hit*> deltaE_contained;
    for (auto front_hit : telescope_candidates.first) {
      if (std::find(telescope_successes.first.begin(),
                    telescope_successes.first.end(),
                    front_hit) == telescope_successes.first.end()) {
        treatDSSSDHit(front_hit);
        deltaE_contained.emplace(front_hit);
        proton_events++;
      }
    }

//    if (proton_events > 2) cout << proton_events << "\t" << endl;
    //Now if the hit in question has EXACTLY 2 proton events we save it as a twoproton event.
    if (proton_events != 2) return;

    vector<Hit*> twoPHits;
    for (auto hit : telescope_successes.first) {
      twoPHits.emplace_back(hit);
      }
    for (auto hit : deltaE_contained) {
      twoPHits.emplace_back(hit);
    }



    double E1 = twoPHits[0]->E;
    double E2 = twoPHits[1]->E;
    Theta = twoPHits[0]->direction.Angle(twoPHits[1]->direction);
    Q2p = E1 + E2 + PROTON_MASS*(E1 + E2 + 2*sqrt(E1*E2)*cos(Theta))/twoPdaughter.getMass();
    Theta *= TMath::RadToDeg();
    // Omega = interp.Eval(Theta);
    Omega = NAN;


    for (auto hit : twoPHits) {
      addTwoProtonHit(hit);
    }
  }

private:
  gsl_vector* theta;
  gsl_vector* omega;
  Interpolator interp = Interpolator(0, kAKIMA);
};




int main(int argc, char *argv[]) {
  clock_t start = clock();
//preparefileio sets up the location of the setup file, target file, input file and output dir
// getbasename takes the name from the config file which we give as command line input followed by run number
cout << "Config file path: " << getProjectRoot() + "analysis/genanalysis/" + getBasename(argv[1]) << endl;
prepareFileIO(getProjectRoot() + "analysis/genanalysis/" + getBasename(argv[1]));

auto setup = JSON::readSetupFromJSON(setup_path);
auto target = make_shared<Target>(JSON::readTargetFromJSON(target_path));

vector<string> input;
int run;

run = stoi(argv[2]); //converts the second command line argument into an integer which is the run number
prepareAnalysis(run); // extracts the specific analysis from the config and the bool values specified in the config

findFilesMatchingWildcard(Form(input_path.c_str(), isotopetype.c_str(),run), input);
      //finds all the run files for a single run number and puts them in input
system(("mkdir -p " + output_path_dir).c_str());




for(auto &runpart : input){
  SortedReader reader{*setup};
  reader.add(runpart);
  reader.setVerbose(true);

  //stem returns "file" of input path/to/file.xyz
  string stem = getStem(runpart);
  TString outfile = (output_path_dir + "/" + stem + "lio.root").c_str(); //specifies the name of the output
  TFile output(outfile, "RECREATE");
  shared_ptr<GeneralAnalysis> analysis;
  shared_ptr<U1analysis> U1ana;
  if(specificAnalysis == "BananaMaker"){
    analysis = make_shared<BananaMaker>(setup, target, &output, isotopetype, exclude_hpges, include_DSSSD_rim,
                                          include_spurious_zone, include_banana_cuts, include_beta_region);
    }//type of analysis can add more else 
  else if(specificAnalysis == "SingleProton"){
    analysis = make_shared<SingleProton>(setup, target, &output, isotopetype, exclude_hpges, include_DSSSD_rim,
                                          include_spurious_zone, include_banana_cuts, include_beta_region);
    }
  else if(specificAnalysis == "GammaSpec"){
    analysis = make_shared<GammaSpec>(setup, target, &output, isotopetype, exclude_hpges, include_DSSSD_rim,
                                          include_spurious_zone, include_banana_cuts, include_beta_region);
    }
  else if(specificAnalysis == "TwoProton"){
    analysis = make_shared<TwoProton>(setup, target, &output, isotopetype, exclude_hpges, include_DSSSD_rim,
                                          include_spurious_zone, include_banana_cuts, include_beta_region);
    }
  else if(specificAnalysis == "AboveBananaAnalysis"){
    U1ana = make_shared<AboveBananaAnalysis>(setup, target, &output, isotopetype, Only_U1);
  }  
  else if(specificAnalysis == "Bananaexplorer"){
    U1ana = make_shared<Bananaexplorer>(setup, target, &output, isotopetype,Only_U1, peaks);
  }  
  else {
      cerr << "Type of analysis not recognizable -- recheck config file -- Aborting analysis." << endl;
      abort();
    }

  printConfig();

  cout << "Reading input from: " << runpart << endl;
  cout << "Writing output to : " << outfile << endl;

  if (analysis) {  // Checks if analysis is not nullptr
    reader.attach(analysis);
  } 
  else if (U1ana) {  // Checks if U1ana is not nullptr
    reader.attach(U1ana);
  }   
  else {
    cout << "No analysis seems to be correctly initialized..." << endl;
  }
  reader.run();
  clock_t stop = clock();
  double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
  printf("\nFile: %s.\tTime elapsed: %.5f\n", (stem + "lio.root").c_str(), elapsed);

  }// end for forloop over runparts
  return 0; // SUCCESS!
  };

