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
string specificAnalysis;
bool exclude_hpges, exclude_U5, include_DSSSD_rim, include_spurious_zone, bananas, exclude_beta_region;

Config cfg;
string path;

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
      output_path_dir = getProjectRoot() + "data/" + getStem(configfile);}
  else{
      output_path_dir = getProjectRoot() + cfg.lookup("output_path_dir").c_str();}
  };//prepareFileIO

void prepareAnalysis(unsigned int run_number){
  specificAnalysis = cfg.lookup("specificAnalysis").c_str();
  exclude_hpges = cfg.exists("exclude_hpges") && cfg.lookup("exclude_hpges");
  exclude_U5 = cfg.exists("exclude_U5") && cfg.lookup("exclude_U5");
  include_DSSSD_rim = cfg.exists("include_dsd_rim") && cfg.lookup("include_DSSSD_rim");
  exclude_beta_region = cfg.exists("include_beta_region") && cfg.lookup("exclude_beta_region");
  include_spurious_zone = cfg.exists("include_spurious_zone") && cfg.lookup("include_spurious_zone");
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
             << "Exclude U5:            " << exclude_U5                                     << endl
             << "Include DSSSD rims:    " << include_DSSSD_rim                              << endl
             << "Exclude beta region:   " << exclude_beta_region                            << endl
             << "Include spurious zone: " << include_spurious_zone                          << endl
             << "-------------------------------------------------------------------------" << endl << endl;
    }
}
#endif