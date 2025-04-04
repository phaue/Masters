#ifndef ANALYSIS_CONFIG_H
#define ANALYSIS_CONFIG_H


#include <string>
#include <cstdio>
#include <libconfig.h++>
#include <ausa/json/IO.h>
#include "projectutil.h"

using namespace std;
using namespace libconfig;
using namespace EUtil;
using namespace AUSA;


string setup_path, target_path, input_path, output_path_dir;
string specificAnalysis, isotopetype;
bool exclude_hpges, include_DSSSD_rim, include_spurious_zone, include_banana_cuts, include_beta_region, Only_U1;
vector<vector<int>> peaks;

Config cfg;
string path;

vector<vector<int>> parsePeaks(const string& peaksStr) {
  vector<vector<int>> peaksSets;
  istringstream iss(peaksStr);
  string set;

  while (getline(iss, set, ';')) {
      vector<int> peaks;
      istringstream setStream(set);
      string token;
      while (getline(setStream, token, ',')) {
          peaks.push_back(stoi(token));
      }
      peaksSets.push_back(peaks);
  }

  return peaksSets;
}


void prepareFileIO(const string& configfile){
  cfg.readFile(configfile.c_str());
  if (cfg.exists("setup_file")) {
    setup_path = getProjectRoot() + cfg.lookup("setup_file").c_str();
  } else {
    cerr << "'setup_file' not found in the config file!" << endl;
  }
  if (cfg.exists("target_file")) {
    target_path = getProjectRoot() + cfg.lookup("target_file").c_str();
  } else {
    cerr << "'target_file' not found in the config file!" << endl;
  }
  if (cfg.exists("input_file")) {
    input_path = getProjectRoot() + cfg.lookup("input_file").c_str();
  } else {
    cerr << "'input_file' not found in the config file!" << endl;
  }
  if (cfg.exists("auto_output_path_dir") && cfg.lookup("auto_output_path_dir")){
      output_path_dir = getProjectRoot() + "data/" + getStem(configfile) +"/" + cfg.lookup("isotopetype").c_str();}
  else{
      output_path_dir = getProjectRoot() + cfg.lookup("output_path_dir").c_str();}

  };//prepareFileIO

void prepareAnalysis(unsigned int run_number){
  specificAnalysis = cfg.lookup("specificAnalysis").c_str();
  isotopetype = cfg.lookup("isotopetype").c_str();
  exclude_hpges = cfg.exists("exclude_hpges") && cfg.lookup("exclude_hpges");
  include_DSSSD_rim = cfg.exists("include_DSSSD_rim") && cfg.lookup("include_DSSSD_rim");
  include_beta_region = cfg.exists("include_beta_region") && cfg.lookup("include_beta_region");
  include_spurious_zone = cfg.exists("include_spurious_zone") && cfg.lookup("include_spurious_zone");
  include_banana_cuts = cfg.exists("include_banana_cuts") && cfg.lookup("include_banana_cuts");
  Only_U1 = cfg.exists("Only_U1") && cfg.lookup("Only_U1");
  if (cfg.exists("peaks")) {
    peaks = parsePeaks(cfg.lookup("peaks").c_str());}

}//prepareAnalysis

void printConfig() {
    if (cfg.lookup("verbose")) {
        cout << "---------------------------- Configuration ------------------------------" << endl
             << "Specific analysis:     " << specificAnalysis                               << endl
             << "Setup:                 " << setup_path                                     << endl
             << "Target:                " << target_path                                    << endl
             << "Input:                 " << input_path                                     << endl
             << "Output:                " << output_path_dir                                << endl
             << "Exclude HPGes:         " << exclude_hpges                                  << endl
             << "Include DSSSD rims:    " << include_DSSSD_rim                              << endl
             << "Exclude beta region:   " << include_beta_region                            << endl
             << "Include spurious zone: " << include_spurious_zone                          << endl
             << "Include_banana_cuts:   " << include_banana_cuts                            << endl
             << "-------------------------------------------------------------------------" << endl << endl;
    }
}
#endif