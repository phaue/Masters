#ifndef U1ANALYSIS_H
#define U1ANALYSIS_H

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
using namespace AUSA::Parser;

class U1analysis : public AbstractSortedAnalyzer{
    public:
      U1analysis(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output,
                      string _isotope, bool _Only_U1)
          :setupSpecs(_setupSpecs), target(_target), isotope(_isotope), Only_U1(_Only_U1){
        //Constructor for the general analysis script

        cout << "input ion is " << ion << endl;


        if(isotope=="Si"){
            implantation_depth = 17/1e6;
        }
        else if(isotope=="P"){
            implantation_depth = 22.4/1e6;
        }
        else if(isotope=="Al"){
            implantation_depth = 41/1e6;
        }
        else if(isotope=="Mg"){
            implantation_depth=45.9/1e6;
        }
        else{
          cout<<"No specification on the input isotope?? Hello?" << endl;
        }

        origin = target->getCenter() + (target->getThickness() / 2. - implantation_depth) * target->getNormal();
       

        //Define all the detectors in the setup using the Detector.h file

        if (!Only_U1){
        U1 = new Detector_frib(0, "U1", DSSSD, Proton, setupSpecs, 500.); //these can be defined with betacutoffs aswell
        U2 = new Detector_frib(1, "U2", DSSSD, Proton, setupSpecs, 500.);
        U3 = new Detector_frib(2, "U3", DSSSD, Proton, setupSpecs, 500.);
        U4 = new Detector_frib(3, "U4", DSSSD, Proton, setupSpecs, 1000.);
        P1 = new Detector_frib(6, "P1", Pad, Proton, setupSpecs);
        P2 = new Detector_frib(7, "P2", Pad, Proton, setupSpecs);
        P3 = new Detector_frib(8, "P3", Pad, Proton, setupSpecs);
        P4 = new Detector_frib(9, "P4", Pad, Proton, setupSpecs);
        makePartners(U1, P1);
        makePartners(U2, P2);
        makePartners(U3, P3);
        makePartners(U4, P4);
        detectors.insert({U1,P1, U2, P2, U3, P3, U4, P4}); 
        }
        else{
          U1 = new Detector_frib(0, "U1", DSSSD, Proton, setupSpecs, 500.); //these can be defined with betacutoffs aswell
          P1 = new Detector_frib(6, "P1", Pad, Alpha, setupSpecs);
          makePartners(U1, P1);
          U1->setBananaCut(new gCut(getProjectRoot() + "data/cuts/totcuts.root", "abovebanU1", include_region));
          detectors.insert({U1,P1}); 

        }
        //Partner the dsssds and pads into telescopes

      
        output->cd(); //something about output file used for mid-analysis dumping
        tree = new TTree("a", "a");
        tree->Branch("num", &NUM);
        tree->Branch("mul", &mul);

        v_id = make_unique<DynamicBranchVector<unsigned short>>(*tree, "id", "mul");

        v_pos = make_unique<DynamicBranchVector<TVector3>>(*tree, "pos");
        v_tarpos = make_unique<DynamicBranchVector<TVector3>>(*tree, "tarpos");
        v_dir = make_unique<DynamicBranchVector<TVector3>>(*tree, "dir");


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

    pAlCalc = defaultRangeInverter("p", "Aluminum");
    pSiCalc = defaultRangeInverter("p", "Silicon"); //Eloss in detector material of protons
    for (auto &layer: target->getLayers()) {
      pTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
      //eloss of protons in target material
    }


    }//General analysis class

    virtual void specificAnalysis() =0;

    //code that runs before we start iterating over events
      void setup(const SortedSetupOutput &output) override{
      AbstractSortedAnalyzer::setup(output); //calls the setup function from the base class NonSpecficAnalysis
      //clock = output.getScalerOutput("clock");
    }//specificanalysis


    void analyze() override {
        clear();
        findHits();
        specificAnalysis();
        tree->Fill();
        NUM++; //counts number of events
      }//analyze

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
          hit.FT = fTime(out, i); // from AUSA
          hit.BT = bTime(out, i); // from AUSA

          TVector3 pos = d.getUniformPixelPosition(FI, BI);
          hit.position = pos;
          TVector3 dir = d.getNormal();
          hit.direction = dir;
          hit.targetposition = NAN_TVECTOR3;
          auto incidenceangle = 0;
          hit.angle = incidenceangle;
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

          auto FI = out.segment(i);
          hit.FI = short(FI);
          hit.FE = hit.Edep;
          hit.FT = out.time(i);

          TVector3 pos = d.getPosition(FI);
          hit.position = pos;
          TVector3 dir = (hit.position).Unit();
          hit.direction = dir;
          auto incidenceangle = 0;
          hit.angle = incidenceangle;

          hits.emplace_back(std::move(hit));
          }//forloop
        }//findPadHit
        bool treatTelescopeHit(Hit *dsssd_hit, Hit *pad_hit) {
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
          double E = pad_hit->Edep;
          if(front_det->getCalibration() == Proton && back_det->getCalibration() == Alpha) {
            E *= 1.016; // This value is extracted from Eriks thesis as a multiplier to somewhat account
                //for the poorer calibrations of the pad detectors.
          }//if statement
          //Energy correction for protons in the
          E+= pSiCalc -> getTotalEnergyCorrection(E, back_det_fdl);
          E+= pAlCalc -> getTotalEnergyCorrection(E, back_det_fct);
          E+= pAlCalc -> getTotalEnergyCorrection(E, front_det_bct);
          E+= pSiCalc -> getTotalEnergyCorrection(E, front_det_bdl);
          E+= dsssd_hit->Edep;
          E+= pSiCalc -> getTotalEnergyCorrection(E, front_det_fdl);
          auto &from = dsssd_hit->position;
          for (auto &intersection: target->getIntersections(from, origin)) {
            auto &calc = pTargetCalcs[intersection.index];
            E += calc->getTotalEnergyCorrection(E, intersection.transversed);
          }//forloop for E energy correction
          //correction with regards to spurious zones
          dsssd_hit->E = E;
          return true;
          }//treatTelescopeHit
    
          bool specialTelescopeTreatment(Hit* dsssd_hit, Hit* pad_hit, double initial_E){
              auto front_det = dsssd_hit->detector;
              auto back_det = pad_hit->detector;

              double angle = 0.0;
              double corrected_E = initial_E;
              auto &from = dsssd_hit->position;
              auto &calc = *pSiCalc;

              for (int iter = 0; iter < 5; ++iter) {
                  double ca = cos(angle);
                  if (abs(ca) < 1e-6) ca = (ca >= 0 ? 1e-6 : -1e-6);
                  double eff_fdl = front_det->getFrontDeadLayer() / abs(ca);
                  
                  for (auto &intersection: target->getIntersections(origin, from)) {
                      auto &tc = pTargetCalcs[intersection.index];
                      corrected_E -= tc->getTotalEnergyLoss(corrected_E, intersection.transversed);
                  }
                  
                  corrected_E = initial_E - pSiCalc->getTotalEnergyLoss(initial_E, eff_fdl);

                  auto func = [&](double* x, double* par) {
                      return TMath::Abs(calc.getTotalEnergyLoss(corrected_E, x[0]) - dsssd_hit->Edep);
                  };
                  calc.getTotalEnergyLoss(corrected_E, 1e10, range);
                  TF1 j("m", func, 0, range, 0);
                  double thickness = j.GetMinimumX();

                  double traversed_thickness = front_det->getThickness()
                                            - front_det->getFrontDeadLayer()
                                            - front_det->getBackDeadLayer()
                                            - front_det->getBackContactThickness();

                  angle = acos(traversed_thickness / thickness);
              }

              auto front_det_fdl = front_det->getFrontDeadLayer()/abs(cos(angle));
              auto front_det_bdl = front_det->getBackDeadLayer()/abs(cos(angle));
              auto front_det_bct = front_det->getBackContactThickness()/abs(cos(angle));
              auto back_det_fct = back_det->getFrontContactThickness()/abs(cos(angle));
              auto back_det_dl = back_det->getFrontDeadLayer()/abs(cos(angle)); 
              double E = pad_hit->Edep;

              E+= pSiCalc -> getTotalEnergyCorrection(E, back_det_dl);
              E+= pAlCalc -> getTotalEnergyCorrection(E, back_det_fct);
              E+= pAlCalc -> getTotalEnergyCorrection(E, front_det_bct);
              E+= pSiCalc -> getTotalEnergyCorrection(E, front_det_bdl);
              E+= dsssd_hit->Edep;
              E+= pSiCalc -> getTotalEnergyCorrection(E, front_det_fdl);
              for (auto &intersection: target->getIntersections(from, origin)) {
                  auto &calc = pTargetCalcs[intersection.index];
                  E += calc->getTotalEnergyCorrection(E, intersection.transversed);
              }
              dsssd_hit->angle = angle;
              dsssd_hit->E = E;
              if(isnan(angle)){
                  double somev = (origin.Z()-dsssd_hit->position.Z())/dsssd_hit->direction.Z();
                  double_t tarx = dsssd_hit->position.X()+somev*dsssd_hit->direction.X();
                  TVector3 tarpos(tarx,dsssd_hit->position.Y(),origin.Z());
                  dsssd_hit->targetposition = tarpos;
              }

              return true;
          }
        void addTelescopeHit(Hit *dsssd_hit, Hit *pad_hit){
            v_id->add(dsssd_hit->id); //id of the dsssd determines the id of the telescope hit
            v_pos->add(dsssd_hit->position);
            v_dir->add(dsssd_hit->direction);
            v_tarpos->add(dsssd_hit->targetposition);
            v_angle->add(dsssd_hit->angle);
            v_Edep->add(NAN);
            v_fEdep->add(dsssd_hit->Edep);
            v_bEdep->add(pad_hit->Edep);
            v_FI->add(dsssd_hit->FI); //front strip hit of the dsssd determines FI
            v_BI->add(dsssd_hit->BI); //back strip hit of the dsssd determines BI
            v_FT->add(dsssd_hit->FT);
            v_BT->add(dsssd_hit->BT);

            v_E->add(dsssd_hit->E);

            mul++;
            }//addTelescopeHit
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


        //code that runs after all events have been analyzed
      void terminate() override { //after analyzing we need to terminate the program
      AbstractSortedAnalyzer::terminate();
      gDirectory->WriteTObject(tree); //writes the output of the analysis to the tree in a file
      }//terminate

      void clear() {
        //We need to clear the memory - done manually in cpp
        //must clear all assigned variables used in analysis before going again
        mul = 0;
        Eg1, Eg2 = NAN;
        AUSA::clear(
        *v_id,
        *v_pos, *v_dir, *v_tarpos,
        *v_angle,
        *v_Edep, *v_fEdep, *v_bEdep,
        *v_FI, *v_BI, *v_FE, *v_BE, *v_FT, *v_BT,
        *v_E
            );
        hits.clear();
        }//clear







ListOfMaterials lom;
Material mat = lom.getMaterial("Silicon"); // get Silicon as default
Ion ion = Ion("H1");
std::shared_ptr<EnergyLossCalculator> calc = defaultRangeInverter(ion, mat); // Range inversion approach (best)
double range=0;
TVector3 origin; // origin of the projectiles to be analyzed
//parameters from the constructor
shared_ptr<Setup> setupSpecs;
shared_ptr<Target> target;
double implantation_depth;
//
TelescopeTabulation *pU1P1;
string isotope;
bool Only_U1;

unordered_set<Detector_frib *> detectors;
Detector_frib *U1, *P1, *U2, *P2,*U3, *P3, *U4, *P4;

TTree *tree; //defines the tree in which we save the new parameters
int NUM;

UInt_t mul{}, CLOCK{}; //TPATTERN{}, TPROTONS{},
SortedSignal clock; //tpattern, tprotons, are these tprotons the time related to the measurements?

Double_t Eg1, Eg2;
unique_ptr<DynamicBranchVector<unsigned short>> v_id;
unique_ptr<DynamicBranchVector<TVector3>> v_pos, v_dir, v_tarpos;
unique_ptr<DynamicBranchVector<double>> v_angle;
unique_ptr<DynamicBranchVector<double>> v_Edep, v_fEdep, v_bEdep;
unique_ptr<DynamicBranchVector<unsigned short>> v_FI, v_BI;
unique_ptr<DynamicBranchVector<double>> v_FE, v_BE;
unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;
unique_ptr<DynamicBranchVector<double>> v_E;



vector<Hit> hits;
unique_ptr<EnergyLossRangeInverter> pSiCalc, pAlCalc;
vector<unique_ptr<EnergyLossRangeInverter>> pTargetCalcs;

}; //Class GeneralAnalysis
#endif