#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>  // For gSystem->GetFromPipe()
#include <iostream>
#include <sstream>
#include "projectutil.h"

using namespace std;
using namespace EUtil;

//specifying input and output directories
string output_dir = getProjectRoot() + "analysis/filelink";


int main(int argc, char *argv[]) {

  string input_dir = getProjectRoot() + argv[1];
  auto *c = new TChain("a");
  cout << "Reading from input directory: " << input_dir << endl;

  string command = "ls " + input_dir + "*.root";
  TString file_list_TString = gSystem->GetFromPipe(command.c_str());
  string file_list = file_list_TString.Data();

  istringstream stream(file_list);
  string file;
  while (stream >> file) {
      string fullPath = input_dir + file;
      cout << "Adding file: " << fullPath << endl;
      c->Add(file.c_str());
  }


  TFile *outputFile = new TFile((output_dir + "/combined.root").c_str(), "RECREATE");
  c->Write();
  outputFile->Close();
  cout << "Chained data written to: " << "combined.root" << endl;

  return 0;
  }