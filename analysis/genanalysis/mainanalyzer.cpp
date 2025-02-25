#include <libconfig.h++>
#include "projectutil.h"
#include "Hit.h"
#include "AnalysisConfig.h"
#include "GeneralAnalysis.h"
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
  BananaMaker(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
                    bool _exclude_hpges = false, bool _exclude_U5 = false,
                    bool _include_DSSSD_rim = false, bool _include_spurious_zone = false,
                    bool _include_banana_cuts=false, bool _include_beta_region =false)
        : GeneralAnalysis(_setupSpecs, _target, output, _exclude_hpges, _exclude_U5,
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
findFilesMatchingWildcard(Form(input_path.c_str(), run), input);
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
  if(specificAnalysis == "BananaMaker"){
    analysis= make_shared<BananaMaker>(setup, target, &output, exclude_hpges, exclude_U5, include_DSSSD_rim,
                                        include_spurious_zone, include_banana_cuts, include_beta_region);
    }//type of analysis can add more else ifs
  else {
      cerr << "Type of analysis not recognizable -- recheck config file -- Aborting analysis." << endl;
      abort();
    }

  printConfig();

  cout << "Reading input from: " << runpart << endl;
  cout << "Writing output to : " << outfile << endl;

  reader.attach(analysis);
  reader.run();
  clock_t stop = clock();
  double elapsed = (double) (stop - start) / CLOCKS_PER_SEC;
  printf("\nFile: %s.\tTime elapsed: %.5f\n", (stem + "lio.root").c_str(), elapsed);

  }// end for forloop over runparts
  return 0; // SUCCESS!
  };

