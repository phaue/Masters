//
// Created by erik on 11/12/23.
//

#include <cxxopts.hpp>
#include <TChain.h>
#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <filesystem>
#include <ausa/json/IO.h>
#include <ausa/eloss/Default.h>
#include <ausa/calibration/DSSDElossCorrector.h>
#include <ausa/eloss/EnergyLossCalculator.h>
#include <telescope/DetectorTabulation.h>

using namespace EUtil;
using namespace std;
using namespace AUSA;
using namespace AUSA::Calibration;
using namespace AUSA::EnergyLoss;
using namespace Telescope;

string output_dir = getProjectRoot() + "/data/cal";

// ./printEnergyCorrections -s ../../setup/setup.json -t ../../setup/target.json -d U1,U2,U3,U6 -E 385.8,905.4,1847,2077,2218 -i 34.6 -y -2.0 > ../../calibration/correctedEnergies.dat
// ./printEnergyCorrections -s ../../setup/setup.json --no-target -I a -d U1,U2 --rim -E 5155.4,5485.74,5804.96 -x -5. -y 0. -z 0.  > ../../calibration/corrected3aEnergiesU1U2.dat
// ./printEnergyCorrections -s ../../setup/setup.json --no-target -I a -d U3,U4 --rim -E 5155.4,5485.74,5804.96 -x 5. -y 0. -z 0.  > ../../calibration/corrected3aEnergiesU3U4.dat
// ./printEnergyCorrections -s ../../setup/setup.json --no-target -I a -d U5 --rim -E 5155.4,5485.74,5804.96 -x 0. -y 5. -z 0.  > ../../calibration/corrected3aEnergiesU5.dat
// ./printEnergyCorrections -s ../../setup/setup.json --no-target -I a -d U6 --rim -E 5155.4,5485.74,5804.96 -x 0. -y -5. -z 0.  > ../../calibration/corrected3aEnergiesU6.dat
int main(int argc, char* argv[]) {
  system(("mkdir -p " + output_dir).c_str());

  /*
   * Define command line options and their default values.
   */
  string setup_path = "setup.json", target_path = "target.json", ion = "p";
  bool include_dssd_rim = false, no_target = false;
  double x = 0., y = 0., z = 0.;
  cxxopts::Options options("printEnergyCorrections",
                           "Calculate and print energy corrections in all detector strips using AUSAlib setup and target files");
  options.add_options()
      ("s,setup", "Path to setup file", cxxopts::value<string>()->default_value(setup_path))
      ("t,target", "Path to target file", cxxopts::value<string>()->default_value(target_path))
      ("d,detector", "Detectors - use ',' as separator, e.g. '-d U1,U2,U3'", cxxopts::value<vector<string>>())
      ("E,energy", "Particle energies (keV) to be corrected - use ',' as separator, e.g. '-E 385.8,905.4,1847'",
          cxxopts::value<vector<double>>())
      ("I,ion", "Ion traversing detector telescope", cxxopts::value<string>()->default_value(ion))
      ("i,implantation", "Ion implantation depth (nm) - Note: Mutually exclusive with z.",
       cxxopts::value<string>()->default_value("Half of target thickness"))
      ("x", "Offset in x direction for source point (mm)", cxxopts::value<double>())
      ("y", "Offset in y direction for source point (mm)", cxxopts::value<double>())
      ("z", "Offset in z direction for source point (mm) - Note: Mutually exclusive with implantation.", cxxopts::value<double>())
      ("r,rim", "Include outermost strips (the rim) on DSSDs (default: false)",
          cxxopts::value<bool>()->default_value("false"))
      ("no-target", "No target (e.g. for 3a source) (default: false) - Note: User must then specify all of x, y and z.",
       cxxopts::value<bool>()->default_value("false"))
       ("h,help", "Print usage")
      ;
  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    cout << options.help() << endl;
    exit(0);
  }
  if (!result.count("detector")) {
    cout << "Must specificy at least one detector." << endl
         << "Try running " << argv[0] << " -h" << endl;
    exit(-1);
  }
  if (!result.count("energy")) {
    cout << "Must specificy at least one energy." << endl
         << "Try running " << argv[0] << " -h" << endl;
    exit(-2);
  }
  if (result.count("no-target") && (!result.count("x") || !result.count("y") || !result.count("z"))) {
    cout << "If there is no target, all of x, y and z must be specified." << endl
         << "Try running " << argv[0] << " -h" << endl;
    exit(-3);
  }
  if (result.count("implantation") && result.count("z")) {
    cout << "implantation and z are mutually exclusive." << endl
         << "Try running " << argv[0] << " -h" << endl;
    exit(-4);
  }

  /*
   * Parse options given on command line.
   */
  include_dssd_rim = result.count("rim") ? result["rim"].as<bool>() : include_dssd_rim;
  no_target = result.count("no-target") ? result["no-target"].as<bool>() : no_target;
  x = result.count("x") ? result["x"].as<double>() : x;
  y = result.count("y") ? result["y"].as<double>() : y;
  z = result.count("y") ? result["y"].as<double>() : y;
  setup_path = result.count("setup") ? result["setup"].as<string>() : setup_path;
  auto setup = JSON::readSetupFromJSON(setup_path);
  target_path = result.count("target") ? result["target"].as<string>() : target_path;
  shared_ptr<Target> target;
  double implantation_depth = 0.;
  TVector3 sourcePos;
  if (!no_target) {
    target = make_shared<Target>(JSON::readTargetFromJSON(target_path));
    implantation_depth = result.count("implantation") ?
                         stod(result["implantation"].as<string>()) : 1e6 * target->getThickness() / 2.; // nm
    implantation_depth /= 1e6; // mm
    sourcePos = target->getCenter()
        + TVector3{x, y, -(target->getThickness()/2. - implantation_depth)};
  } else {
    target_path = "NOT USING TARGET";
    sourcePos = TVector3{x, y, z};
  }
  ion = result.count("ion") ? result["ion"].as<string>() : ion;


  /*
   * Logging.
   */
  string pwd = std::filesystem::current_path();
  time_t now = time(nullptr);
  char* datetime = ctime(&now);
  cout << "# " << datetime
       << "# Output created from within " << pwd << " with the following command" << endl << "# ";
  for (int i = 0; i < argc; i++) {
    cout << argv[i] << " ";
  }
  cout << endl
       << "# Setup:          " << setup_path << endl
       << "# Target:         " << target_path << endl
       << "# Source point:   " << "(" << sourcePos.X() << ", " << sourcePos.Y() << ", "
                                      << sourcePos.Z() << ")"  << endl
       << "# Energies (keV): ";
  for (const auto& E : result["energy"].as<vector<double>>()) {
    cout << E << " ";
  }
  cout << endl;

  for (const auto& det : result["detector"].as<vector<string>>()) {
    auto dT = new DetectorTabulation(setup, target, det, ion);
    if (!no_target) dT->setImplantationDepth(implantation_depth);
    
    auto fCount = dT->getFrontStripCount();
    auto bCount = dT->getBackStripCount();
    auto dl = dT->getFrontDeadLayerThickness();

    auto corrector = dT->getDetectorEnergyLossCorrector();

    cout << "# Detector=" << det << "\t"
         << "FRONT STRIPS " << 1 + !include_dssd_rim << ".." << fCount - !include_dssd_rim << endl;
    
    vector<UInt_t> disabled;
    if (!include_dssd_rim) { disabled.push_back(1); disabled.push_back(16);  }
    for (int i = 0 + !include_dssd_rim; i < fCount - !include_dssd_rim; i++) {
      for (const auto& E : result["energy"].as<vector<double>>()) {
        cout << dT->getWeightedStripDepositedEnergy(E, front, i, disabled) << "\t";
      }
      cout << "# " << i+1 << endl;
    }
    cout << "# Detector=" << det << "\t"
         << "BACK STRIPS " << 1 + !include_dssd_rim << ".." << bCount - !include_dssd_rim << endl;
    for (int i = 0 + !include_dssd_rim; i < bCount - !include_dssd_rim; i++) {
      for (const auto& E : result["energy"].as<vector<double>>()) {
        cout << dT->getWeightedStripDepositedEnergy(E, back, i, disabled) << "\t";
      }
      cout << "# " << i+1 << endl;
    }
    cout << endl;
  }

  return 0;
}