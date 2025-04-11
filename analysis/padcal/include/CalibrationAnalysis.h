#ifndef GENERALANALYSIS_H
#define GENERALANALYSIS_H

#include <libconfig.h++>
#include "projectutil.h"
#include "Hit.h"
#include "Detector_frib.h"

#include <TROOT.h>

#include <ausa/json/IO.h>
#include <ausa/setup/Target.h>
#include <ausa/eloss/Default.h>
#include <ausa/util/FileUtil.h>
#include <ausa/sort/SortedReader.h>
#include <ausa/util/DynamicBranchVector.h>
#include <ausa/output/OutputConvenience.h>
#include <ausa/sort/analyzer/AbstractSortedAnalyzer.h>
#include <telescope/TelescopeTabulation.h>
#include <ausa/constants/Mass.h>
#include <ausa/eloss/Ion.h>

#include <ausa/json/MaterialParser.h>
#include <ausa/json/TargetParser.h>
#include <ausa/parser/UnitParser.h>
#include <ausa/parser/ThicknessParser.h>
#include <ausa/eloss/Material.h>
#include <ausa/eloss/ListOfMaterials.h>
#include <ausa/eloss/EnergyLossCalculator.h>
#include <ausa/eloss/Default.h>

#include <TVector3.h>

#include <unistd.h>
#include <getopt.h>
#include <TF1.h>

#define NAN_UINT 100
#define NAN_TVECTOR3 TVector3(NAN, NAN, NAN)

using namespace std;
using namespace libconfig;
using namespace EUtil;
using namespace Detectors_;

using namespace AUSA;
using namespace AUSA::Sort;
using namespace AUSA::EnergyLoss;
using namespace Telescope;
using namespace AUSA::Constants;

class CalibrationAnalysis : public AbstractSortedAnalyzer{
    public:
      CalibrationAnalysis(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output, string _isotopetype, bool _events_matched=false)
          :setupSpecs(_setupSpecs), target(_target), outFile(output), isotopetype(_isotopetype), events_matched(_events_matched){
        //Constructor for the general analysis script

        if(isotopetype=="Si"){
          implantation_depth = 17/1e6;
      }
      else if(isotopetype=="P"){
          implantation_depth = 22.4/1e6;
      }
      else if(isotopetype=="Al"){
          implantation_depth = 41/1e6;
      }
      else if(isotopetype=="Mg"){
          implantation_depth=45.9/1e6;
      }
      else{
        cout<<"No specification on the input isotope?? Hello?" << endl;
      }

      origin = target->getCenter() + (target->getThickness() / 2. - implantation_depth) * target->getNormal();

        U1 = new Detector_frib(0, "U1", DSSSD, Proton, setupSpecs, 500.); //these can be defined with betacutoffs aswell
        U2 = new Detector_frib(1, "U2", DSSSD, Proton, setupSpecs, 500.);
        U3 = new Detector_frib(2, "U3", DSSSD, Proton, setupSpecs, 500.);
        P1 = new Detector_frib(6, "P1", Pad, NoCalibration, setupSpecs);
        P2 = new Detector_frib(7, "P2", Pad, NoCalibration, setupSpecs);
        P3 = new Detector_frib(8, "P3", Pad, NoCalibration, setupSpecs);
        makePartners(U1, P1);
        makePartners(U2, P2);
        makePartners(U3, P3);
        detectors.insert({U1, U2, U3, P1, P2, P3}); 

        outFile->cd(); //something about output file used for mid-analysis dumping
        tree = new TTree("a", "a");
        tree->Branch("num", &NUM);
        tree->Branch("mul", &mul);

        v_id = make_unique<DynamicBranchVector<unsigned short>>(*tree, "id", "mul");

        v_dir = make_unique<DynamicBranchVector<TVector3>>(*tree, "dir");
        v_pos = make_unique<DynamicBranchVector<TVector3>>(*tree, "pos");

        v_theta = make_unique<DynamicBranchVector<double>>(*tree, "theta", "mul");
        v_phi = make_unique<DynamicBranchVector<double>>(*tree, "phi", "mul");
        v_angle = make_unique<DynamicBranchVector<double>>(*tree, "angle", "mul"); // angle of incidence w.r.t. detector surface

        v_Edep = make_unique<DynamicBranchVector<double>>(*tree, "Edep", "mul");
        v_fEdep = make_unique<DynamicBranchVector<double>>(*tree, "fEdep", "mul");
        v_bEdep = make_unique<DynamicBranchVector<double>>(*tree, "bEdep", "mul");

        v_FI = make_unique<DynamicBranchVector<unsigned short>>(*tree, "FI", "mul");
        v_BI = make_unique<DynamicBranchVector<unsigned short>>(*tree, "BI", "mul");
        v_FE = make_unique<DynamicBranchVector<double>>(*tree, "FE", "mul");
        v_BE = make_unique<DynamicBranchVector<double>>(*tree, "BE", "mul");
        v_FT = make_unique<DynamicBranchVector<double>>(*tree, "FT", "mul");
        v_BT = make_unique<DynamicBranchVector<double>>(*tree, "BT", "mul");

        v_E = make_unique<DynamicBranchVector<double>>(*tree, "E", "mul");
        v_Ech = make_unique<DynamicBranchVector<double>>(*tree, "Ech", "mul");
        v_Ecal = make_unique<DynamicBranchVector<double>>(*tree, "Ecal", "mul");


        if (!events_matched){
        histP1 = make_unique<TH1D>("P1ch", "P1ch", 2500,0,2500);
        histP2 = make_unique<TH1D>("P2ch", "P2ch", 2500,0,2500);
        histP3 = make_unique<TH1D>("P3ch", "P3ch", 2500,0,2500);
        }
        else {
        histP1 = make_unique<TH1D>("P1cal", "P1cal", 700,0,7000);
        histP2 = make_unique<TH1D>("P2cal", "P2cal", 700,0,7000);
        histP3 = make_unique<TH1D>("P3cal", "P3cal", 700,0,7000);
        }



    pSiCalc = defaultRangeInverter("p", "Silicon"); //Eloss in detector material of protons
    for (auto &layer: target->getLayers()) {
      pTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
      //eloss of protons in target material
    }

    pAlCalc = defaultRangeInverter("p", "Aluminum"); //Eloss in detector material of protons
    for (auto &layer: target->getLayers()) {
      pTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
      //eloss of protons in target material
    }
    }//General analysis class constructors


    virtual void specificAnalysis() =0;

    //code that runs before we start iterating over events
      void setup(const SortedSetupOutput &output) override{
      AbstractSortedAnalyzer::setup(output); 
    }//specificanalysis


    void analyze() override {
        clear();
        findHits();
        specificAnalysis();
        tree->Fill();
        NUM++;
      }//analyze

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Find hits & Treat hits functions are put here since they are related to the analyze function, although not directly



    void findHits() {
        for(const auto &det: detectors) {
          auto type = det->getType();
          switch (type) {
            case DSSSD:
              findDSSSDHit(det);
              break;
            case Pad:
              findPadHit(det);
              break;
             case HPGe:
              break;
             case NoType:
               cerr << "!WARNING! --> The detector " << det->getName() << " is of type NoType <-- !WARNING!"  << endl;
               break;
              }//switch
          }//forloop
      }//findHits


      void findDSSSDHit(Detector_frib *detector) {
        unsigned short id = detector->getId();
        auto &out = output.getDssdOutput(detector->getName());
        auto &d = out.detector();
        auto MUL = AUSA::mul(out);

        for (int i = 0; i<MUL; i++) {
          Hit hit;
          hit.detector = detector;
          hit.id = id;
          hit.Edep = energy(out, i); // from AUSA


          auto FI = fSeg(out, i); // from AUSA
          auto BI = bSeg(out, i); // from AUSA
        
          hit.FI = short(FI);
          hit.BI = short(BI);
          hit.FE = fEnergy(out, i); // from AUSA
          hit.BE = bEnergy(out, i); // from AUSA

          hits.emplace_back(std::move(hit));
          }//forloop
        }//find dsssd hit


      void findPadHit(Detector_frib *detector) {
        unsigned short id = detector->getId();
        auto &out = output.getSingleOutput(detector->getName());
        auto &d = out.detector();
        auto MUL = AUSA::mul(out);

        for (int i = 0; i<MUL; i++) {
          Hit hit;
          hit.detector = detector;
          hit.id = id;
          hit.Edep = out.energy(i); // from AUSA
          hit.Ech = hit.Edep;
          auto FI = out.segment(i);
          hit.FI = short(FI);
          hit.FE = hit.Edep;

          hits.emplace_back(std::move(hit));
          }//forloop
        }//findPadHit


        bool treatBackwardHit(Hit *dsssd_hit, Hit *pad_hit, double peak_E){
          auto front_det = dsssd_hit->detector;
          auto back_det = pad_hit->detector;
          double angle = dsssd_hit->angle;
          auto front_det_fdl = front_det->getFrontDeadLayer()/abs(cos(angle));
          auto front_det_bdl = front_det->getBackDeadLayer()/abs(cos(angle));
          auto front_det_bct = front_det->getBackContactThickness()/abs(cos(angle));
          auto back_det_fct = back_det->getBackContactThickness()/abs(cos(angle));
          auto back_det_fdl = back_det->getFrontDeadLayer()/abs(cos(angle));

          double pad_E = peak_E;
          auto &from = dsssd_hit->position;
          for (auto &intersection: target->getIntersections(origin, from)) {
            auto &calc = pTargetCalcs[intersection.index];
            pad_E += calc->getTotalEnergyLoss(peak_E, intersection.transversed);
          }//forloop for finding the original particle energy
          pad_E -= pSiCalc -> getTotalEnergyLoss(pad_E, front_det_fdl);
          pad_E -= dsssd_hit->Edep;
          pad_E -= pSiCalc -> getTotalEnergyLoss(pad_E, front_det_bdl);
          pad_E -= pAlCalc -> getTotalEnergyLoss(pad_E, front_det_bct);
          pad_E -= pAlCalc -> getTotalEnergyLoss(pad_E, back_det_fct);
          pad_E -= pSiCalc -> getTotalEnergyLoss(pad_E, back_det_fdl);

          pad_hit->Ecal = pad_E;
          return true;
        }

        bool treatTelescopeHit(Hit *dsssd_hit, Hit *pad_hit) {

          /*
          This does not seem to work, wonder if i am misunderstanding the functions here
          */

          // if there is no energy recorded in either then there is no telescope hit therefore return false
          if (dsssd_hit->Edep == 0 || pad_hit->Edep == 0) return false;

          auto front_det = dsssd_hit->detector;
          auto back_det = pad_hit->detector;
          double angle = dsssd_hit->angle;

          auto front_det_fdl = front_det->getFrontDeadLayer()/abs(cos(angle));
          auto front_det_bdl = front_det->getBackDeadLayer()/abs(cos(angle));
          auto front_det_bct = front_det->getBackContactThickness()/abs(cos(angle));
          auto back_det_fct = back_det->getBackContactThickness()/abs(cos(angle));
          auto back_det_fdl = back_det->getFrontDeadLayer()/abs(cos(angle));
          
          auto active_thickness = (front_det->getThickness()-front_det->getFrontDeadLayer()-front_det->getBackDeadLayer()-front_det->getBackContactThickness())/abs(cos(angle));
          double iniE = dsssd_hit->Edep + pSiCalc -> getTotalEnergyCorrection(dsssd_hit->Edep, active_thickness);
          double origE = iniE + pSiCalc -> getTotalEnergyCorrection(iniE, front_det_fdl);
          /*auto &from = dsssd_hit->position;
          for (auto &intersection: target->getIntersections(from, origin)) {
            auto &calc = pTargetCalcs[intersection.index];
            origE += calc->getTotalEnergyCorrection(origE, intersection.transversed);
          }//forloop for finding the original particle energy

          iniE -= pSiCalc -> getTotalEnergyLoss(iniE, front_det_bdl);
          iniE -= pAlCalc -> getTotalEnergyLoss(iniE, front_det_bct);
          iniE -= pAlCalc -> getTotalEnergyLoss(iniE, back_det_fct);
          iniE -= pSiCalc -> getTotalEnergyLoss(iniE, back_det_fdl);
          */
          
          dsssd_hit->Edep = origE;
          pad_hit->Edep = iniE;
          return true;
          }//treatTelescopeHit
          void addTelescopeHit(Hit *dsssd_hit, Hit *pad_hit){
            v_id->add(dsssd_hit->id); //id of the dsssd determines the id of the telescope hit
            v_dir->add(dsssd_hit->direction);
            v_pos->add(dsssd_hit->position);
            v_theta->add(dsssd_hit->theta);
            v_phi->add(dsssd_hit->phi);
            v_angle->add(dsssd_hit->angle);
            v_Edep->add(NAN);
            v_fEdep->add(dsssd_hit->Edep);
            v_bEdep->add(pad_hit->Edep);
            v_FI->add(dsssd_hit->FI); 
            v_BI->add(dsssd_hit->BI);
            v_Ech->add(pad_hit->Ech);
            v_E->add(dsssd_hit->E);
            v_Ecal->add(pad_hit->Ecal);
            mul++;
            }//addTelescopeHit
            void addPadHit(Hit *hit){
              v_id->add(hit->id);
              v_Edep->add(hit->Edep);
              v_E->add(hit->E);
              mul++;
              //p = true
            }//addPadHit
        //code that runs after all events have been analyzed
      void terminate() override { //after analyzing we need to terminate the program
      AbstractSortedAnalyzer::terminate();
      outFile->cd();
      //Write all the histograms and the TTree object to the output file.
      tree->Write();
      histP1->Write();
      histP2->Write();
      histP3->Write();
      outFile->Close();
      }//terminate

      void clear() {
        mul = 0;
        AUSA::clear(
        *v_id,
        *v_dir, *v_pos,
        *v_theta, *v_phi, *v_angle,
        *v_Edep, *v_fEdep, *v_bEdep,
        *v_FI, *v_BI, *v_FE, *v_BE, *v_FT, *v_BT,
        *v_E, *v_Ech, *v_Ecal
            );
        hits.clear();
        }//clear




bool events_matched;
unique_ptr<TH1D> histP1, histP2, histP3;
shared_ptr<Setup> setupSpecs;
shared_ptr<Target> target;
double implantation_depth;
string isotopetype;
//
unordered_set<Detector_frib *> detectors;
Detector_frib *U1, *U2, *U3, *P1, *P2, *P3;//, *U5, *U6, *P5, *P6

TVector3 origin;
TTree *tree; //defines the tree in which we save the new parameters
int NUM;
UInt_t mul{}, CLOCK{}; //TPATTERN{}, TPROTONS{},
SortedSignal clock; //tpattern, tprotons, are these tprotons the time related to the measurements?

unique_ptr<DynamicBranchVector<unsigned short>> v_id;
unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
unique_ptr<DynamicBranchVector<double>> v_theta, v_phi, v_angle;
unique_ptr<DynamicBranchVector<double>> v_Edep, v_fEdep, v_bEdep;
unique_ptr<DynamicBranchVector<unsigned short>> v_FI, v_BI;
unique_ptr<DynamicBranchVector<double>> v_FE, v_BE;
unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;
unique_ptr<DynamicBranchVector<double>> v_E, v_Ecal, v_Ech;

vector<Hit> hits;
unique_ptr<EnergyLossRangeInverter> pSiCalc;
unique_ptr<EnergyLossRangeInverter> pAlCalc;
vector<unique_ptr<EnergyLossRangeInverter>> pTargetCalcs;
TFile* outFile;

}; //Class GeneralAnalysis
#endif