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

bool events_matched;
string setup_path, target_path, input_path, output_path_dir;
string specificAnalysis, isotopetype;

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
      output_path_dir = getProjectRoot() + "data/" + getStem(configfile) +"/" + cfg.lookup("isotopetype").c_str();}
  else{
      output_path_dir = getProjectRoot() + cfg.lookup("output_path_dir").c_str();}

  };//prepareFileIO

void prepareAnalysis(){
  specificAnalysis = cfg.lookup("specificAnalysis").c_str();
  isotopetype = cfg.lookup("isotopetype").c_str();
  events_matched = cfg.exists("events_matched") && cfg.lookup("events_matched");
}//prepareAnalysis

void printConfig() {
    if (cfg.lookup("verbose")) {
        cout << "---------------------------- Configuration ------------------------------" << endl
             << "Specific analysis:     " << specificAnalysis                               << endl
             << "Isotope type:          " << isotopetype                                    << endl
             << "Setup:                 " << setup_path                                     << endl
             << "Target:                " << target_path                                    << endl
             << "Input:                 " << input_path                                     << endl
             << "Output:                " << output_path_dir                                << endl
             << "Events_matched:        " << events_matched                                 << endl << endl;
            }
}
#endif