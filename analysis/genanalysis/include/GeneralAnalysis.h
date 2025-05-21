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

class GeneralAnalysis : public AbstractSortedAnalyzer{
    public:
      GeneralAnalysis(const shared_ptr<Setup> &_setupSpecs, const shared_ptr<Target> &_target, TFile *output, string _isotopetype,
                      bool _exclude_hpges = false, bool _include_DSSSD_rim = false, bool _include_spurious_zone = false,
                      bool _include_banana_cuts=false, bool _include_beta_region =false)
          :setupSpecs(_setupSpecs), target(_target), isotopetype(_isotopetype),
            exclude_hpges(_exclude_hpges), include_DSSSD_rim(_include_DSSSD_rim), include_spurious_zone(_include_spurious_zone),
             include_banana_cuts(_include_banana_cuts), include_beta_region(_include_beta_region){
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
      //  cout << "target center=(" << target->getCenter().X() << ", " << target->getCenter().Y() << ", " << target->getCenter().Z() << ")" << endl
      //   << " implantation=(" << origin.X() << ", " << origin.Y() << ", " << origin.Z() << ")" << endl
      //   << "target thickness=" << target->getThickness() << endl;

        //Define all the detectors in the setup using the Detector.h file
              //betacutoff can now be determined for U1 U2 U3 U4
        U1 = new Detector_frib(0, "U1", DSSSD, Proton, setupSpecs, 500.); //these can be defined with betacutoffs aswell
        U2 = new Detector_frib(1, "U2", DSSSD, Proton, setupSpecs, 400.);
        U3 = new Detector_frib(2, "U3", DSSSD, Proton, setupSpecs, 350.);
        U4 = new Detector_frib(3, "U4", DSSSD, Proton, setupSpecs, 750.);// should change according to thicknesses
        U5 = new Detector_frib(4, "U5", DSSSD, Alpha, setupSpecs, 1500.);
        U6 = new Detector_frib(5, "U6", DSSSD, Alpha, setupSpecs, 400.);

        P1 = new Detector_frib(6, "P1", Pad, Alpha, setupSpecs);
        P2 = new Detector_frib(7, "P2", Pad, Alpha, setupSpecs);
        P3 = new Detector_frib(8, "P3", Pad, Proton, setupSpecs);
        P4 = new Detector_frib(9, "P4", Pad, Alpha, setupSpecs);
        P5 = new Detector_frib(10, "P5", Pad, Alpha, setupSpecs); 
        P6 = new Detector_frib(11, "P6", Pad, Alpha, setupSpecs);

        G1 = new Detector_frib(12, "G1", HPGe, Gamma, setupSpecs);
        G2 = new Detector_frib(13, "G2", HPGe, Gamma, setupSpecs);

        //Partner the dsssds and pads into telescopes
        makePartners(U1, P1);
        makePartners(U2, P2);
        makePartners(U3, P3);
        makePartners(U4, P4);
        makePartners(U6, P6);

        if (include_banana_cuts) {
        // need to set the banana cuts here when they are done example of how this is done is seen below
          try {
              U1->setBananaCut(new gCut(getProjectRoot() + "data/cuts/banana_cuts.root", "bananaU1", include_region));
              U2->setBananaCut(new gCut(getProjectRoot() + "data/cuts/banana_cuts.root", "bananaU2", include_region));
              U3->setBananaCut(new gCut(getProjectRoot() + "data/cuts/banana_cuts.root", "bananaU3", include_region));
              U4->setBananaCut(new gCut(getProjectRoot() + "data/cuts/banana_cuts.root", "bananaU4", include_region));
              
              U6->setBananaCut(new gCut(getProjectRoot() + "data/cuts/U6cut.root", "CUTG", include_region));

            } catch (const runtime_error &e) {
              cerr << "Error initializing banana cuts: " << e.what() << endl;
              throw;
          }
      } // bananas

        //beta region is the region that lies in the very low energies, these zones are different for each detector
        //i should ideally include them, change detector header and this file when done

        //set telescope tabulations -->> to be seen whether or not these are required remove them if they have no relevance
        pU1P1 = new TelescopeTabulation(setupSpecs, target, "U1", "P1", "p");
        pU2P2 = new TelescopeTabulation(setupSpecs, target, "U2", "P2", "p");
        pU3P3 = new TelescopeTabulation(setupSpecs, target, "U3", "P3", "p");
        pU4P4 = new TelescopeTabulation(setupSpecs, target, "U4", "P4", "p");
        pU6P6 = new TelescopeTabulation(setupSpecs, target, "U6", "P6", "p");


        //Sets the implantation depth for all the tabulations missing , pU6P6, aU6P6
          for (auto &tabulations : {pU1P1, pU2P2, pU3P3, pU4P4, pU6P6}) {
          tabulations->setImplantationDepth(implantation_depth);}
        
        pU1P1->setESignalThreshold(280.);   
        pU2P2->setESignalThreshold(340.); 
        pU3P3->setESignalThreshold(180.); 
        pU4P4->setESignalThreshold(350.); 
        pU6P6->setESignalThreshold(220.);

        U1->addTelescopeTabulation(pU1P1);
        U2->addTelescopeTabulation(pU2P2);
        U3->addTelescopeTabulation(pU3P3);
        U4->addTelescopeTabulation(pU4P4);
        U6->addTelescopeTabulation(pU6P6);

//missing U5, U6, P5, P6
        detectors.insert({U1, U2, U3, U4, U5, U6, P1, P2, P3, P4, P5, P6, G1, G2}); 
        //detectors.insert({U1, U2, U3, U4, P1, P2, P3, P4, G1, G2}); 

        output->cd(); //something about output file used for mid-analysis dumping
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
        v_Ea = make_unique<DynamicBranchVector<double>>(*tree, "Ea", "mul");


        //tree->Branch("Theta", &Theta); // wat the fuck is this
        //tree->Branch("Omega", &Omega);// what the fuck is this

        tree->Branch("Eg1", &Eg1); 
        tree->Branch("Eg2", &Eg2); 

        tree->Branch("pg1", &pg1); 
        tree->Branch("pg2", &pg2); 
        tree->Branch("peak", &peakval);
        tree->Branch("Egated", &Egated);
        //tree->Branch("bg", &bg);
        //tree->Branch("CLOCK", &CLOCK);

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

    aSiCalc = defaultRangeInverter("a", "Silicon");
    for (auto &layer: target->getLayers()) {
      aTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("a"), layer.getMaterial()));
    }


    }//General analysis class


    virtual void specificAnalysis() =0;

    //code that runs before we start iterating over events
      void setup(const SortedSetupOutput &output) override{
      AbstractSortedAnalyzer::setup(output); //calls the setup function from the base class NonSpecficAnalysis
      //clock = output.getScalerOutput("clock");
    }//specificanalysis

    //code that runs for each event in the file this means that each event carries some info about
    //Dsssd hits, pad hits and gamma hits, meaning they are grouped already
    void analyze() override {
      // analyze function, first it clears, then it starts the clock
      //then the hits are found and the analysis run and the tree filled
        clear();
        //CLOCK = clock.getValue();
        findHits();
        specificAnalysis();
        pg1 = p && g1;
        pg2 = p && g2;
        //bg = b && g;
        tree->Fill();
        NUM++; //counts number of events
      }//analyze

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Find hits & Treat hits functions are put here since they are related to the analyze function, although not directly



    void findHits() {
        for(const auto &det: detectors) {
          if (exclude_hpges && det->getType() == HPGe) continue;
          auto type = det->getType();
          switch (type) {
            case DSSSD:
              findDSSSDHit(det);
              break;
            case Pad:
              findPadHit(det);
              break;
             case HPGe:
              findGermaniumHit(det);
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

          if(!include_beta_region && hit.Edep <= detector->getBetaCut()) continue;

          auto FI = fSeg(out, i); // from AUSA
          auto BI = bSeg(out, i); // from AUSA
          if (!include_DSSSD_rim && (FI == 1 || FI == 16 || BI == 1 || BI == 16)) continue;
          if (!include_spurious_zone) {
            double Energy_threshold = INFINITY;
            for (auto &tabulation : detector->getTelescopeTabulations()) {
              if (tabulation->getIon() != "p") continue;
              if (Energy_threshold > tabulation->getEnergyDepositionAtReachThrough(FI, BI)) {
                Energy_threshold = tabulation->getEnergyDepositionAtReachThrough(FI, BI);
              }//if statement
            }//forloop
            // this part here checks whether or not the Edep energy lies within the zone where we cannot "see" the energy
            // due to the energy lying between the reach through energy and punch through energy, if the particle stops in either
            //the back deadlayer of the DSSSD or the front dead layer of the pad.
            if (hit.Edep >= Energy_threshold) continue;
          }//spurious zone correction

          hit.FI = short(FI);
          hit.BI = short(BI);
          hit.FE = fEnergy(out, i); // from AUSA
          hit.BE = bEnergy(out, i); // from AUSA
          hit.FT = fTime(out, i); // from AUSA
          hit.BT = bTime(out, i); // from AUSA

          TVector3 pos = d.getUniformPixelPosition(FI, BI);
          hit.position = pos;
          TVector3 dir = (hit.position - origin).Unit();
          hit.direction = dir;

          hit.theta = dir.Theta();
          hit.phi = dir.Phi();
          auto incidenceangle = hit.direction.Angle(-d.getNormal());
          hit.angle = incidenceangle;
          hits.emplace_back(std::move(hit));
          }//forloop
        }//find dsssd hit


      void findPadHit(Detector_frib *detector) {
        if (detector->getName() == "P5") return; // Pad still dead
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
          TVector3 dir = (hit.position - origin).Unit();
          hit.direction = dir;

          hit.theta = dir.Theta();
          hit.phi = dir.Phi();
          auto incidenceangle = hit.direction.Angle(-d.getNormal());
          hit.angle = incidenceangle;

          hits.emplace_back(std::move(hit));
          }//forloop
        }//findPadHit

        void findGermaniumHit(Detector_frib *detector) {
          unsigned short id = detector->getId();
          auto &out = output.getSingleOutput(detector->getName());
          auto &d = out.detector();
          auto MUL = AUSA::mul(out);
          
          for (int i=0;i<MUL; i++){
            if (id == G1->getId()){
              Eg1 = out.energy(i);
              g1 = true;
            }
            if (id == G2->getId()){
              Eg2 = out.energy(i);
              g2=true;
            }
          }
        }

/*This function essentially sets the energy of the DSSSD hits to its true value by accounting for losses
            in the dead layer of the detector and the target
 */
        void treatDSSSDHit(Hit *hit) {
          auto det = hit->detector;
          double angle = hit->angle;
          //finds the front dead layer accounting for the thickness change with changing angle of incidence
          double fdl = det->getFrontDeadLayer()/abs(cos(angle));

          double E = hit->Edep; // proton energy correction from the front dead layer and target layers
          E += pSiCalc -> getTotalEnergyCorrection(E, fdl);
          auto &from = hit->position;
          for (auto &intersection: target->getIntersections(from, origin)) {
            auto &calc = pTargetCalcs[intersection.index];
            E += calc->getTotalEnergyCorrection(E, intersection.transversed);
          }//forloop for E energy correction
          hit->E = E; // set the energy of the hit to this energy corrected value
          
          double Ea = hit->Edep;
          Ea += aSiCalc->getTotalEnergyCorrection(Ea, fdl);
          for (auto &intersection: target->getIntersections(from, origin)) {
              auto &calc = aTargetCalcs[intersection.index];
              Ea += calc->getTotalEnergyCorrection(Ea, intersection.transversed);
          }
          hit->Ea = Ea;

        }//treatDSSSDHit

        void treatPadHit(Hit *hit){
          auto det = hit->detector;
          double angle = hit->angle;
          double fct = det->getFrontContactThickness()/abs(cos(angle));
          double fdl = det->getFrontDeadLayer()/abs(cos(angle));

          double E = hit->Edep;
          E += pAlCalc -> getTotalEnergyCorrection(E, fct);
          E += pSiCalc -> getTotalEnergyCorrection(E, fdl);
        auto &from = hit->position;
        for (auto &intersection: target->getIntersections(from, origin)) {
          auto &calc = pTargetCalcs[intersection.index];
          E += calc->getTotalEnergyCorrection(E, intersection.transversed);
        }//forloop for E energy correction
        hit->E = E;
        }//treatPadHit

        bool treatTelescopeHit(Hit *dsssd_hit, Hit *pad_hit) {
          // if there is no energy recorded in either then there is no telescope hit therefore return false
          if (dsssd_hit->Edep == 0 || pad_hit->Edep == 0) return false;
          if (dsssd_hit->detector->getName()=="U4"){
          if (dsssd_hit->Edep<dsssd_hit->detector->getBetaCut() && pad_hit->Edep<dsssd_hit->detector->getBetaCut()*1.5) {
            b = true;
          }}

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
              E*= 1.014;
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
          if (!include_spurious_zone) {
            auto FI = dsssd_hit->FI;
            auto BI = dsssd_hit->BI;
//I think these values are to ensure there are no ambiguities with the punch through or reach through values
            double E_lower_threshold = INFINITY, E_upper_threshold = -1.*INFINITY;
            for (auto &tabulation : front_det->getTelescopeTabulations()) {
              if (tabulation->getIon() != "p") continue;
              if(E_lower_threshold > tabulation -> getParticleEnergyAtPunchThrough(FI, BI)){
                E_lower_threshold = tabulation -> getParticleEnergyAtPunchThrough(FI, BI);
                } //if
              if(E_upper_threshold < tabulation -> getParticleEnergyAtReachThrough(FI, BI)){
                E_upper_threshold = tabulation -> getParticleEnergyAtReachThrough(FI, BI);
              } //if
              }//forloop
          if (E_lower_threshold <= E && E <= E_upper_threshold) {
            return false;
          } // this ensures that spurious zones are omitted from the telscope hits, since they cannot theoretically happen
          //since energies between punch through and reach through energies are not able to be seen.
          //is this a problem when determining thicknesses?
          }//if spurious zones
          dsssd_hit->E = E;
          return true;
          }//treatTelescopeHit


//adding a dsssd hit by assigning the corresponding values
          void addDSSSDHit(Hit *hit){
            v_id->add(hit->id);
            v_dir->add(hit->direction);
            v_pos->add(hit->position);
            v_theta->add(hit->theta);
            v_phi->add(hit->phi);
            v_angle->add(hit->angle);
            v_Edep->add(hit->Edep);
            v_fEdep->add(NAN);
            v_bEdep->add(NAN);
            v_FI->add(hit->FI);
            v_BI->add(hit->BI);
            v_FT->add(hit->FT);
            v_BT->add(hit->BT);

            v_E->add(hit->E);
            v_Ea->add(hit->Ea);
            mul++;
            p = true;
            }//addDSSSDHit

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
            v_FI->add(dsssd_hit->FI); //front strip hit of the dsssd determines FI
            v_BI->add(dsssd_hit->BI); //back strip hit of the dsssd determines BI
            v_FT->add(dsssd_hit->FT);
            v_BT->add(dsssd_hit->BT);

            v_E->add(dsssd_hit->E);
            mul++;
            p = true;
            }//addTelescopeHit
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        void addPadHit(Hit *hit){
          v_id->add(hit->id);
          v_dir->add(hit->direction);
          v_pos->add(hit->position);
          v_theta->add(hit->theta);
          v_phi->add(hit->phi);
          v_angle->add(hit->angle);
          v_Edep->add(hit->Edep);
          v_fEdep->add(NAN);
          v_bEdep->add(NAN);
          v_FI->add(hit->FI);
          v_BI->add(hit->BI);
          v_FT->add(hit->FT);
          v_BT->add(hit->BT);

          v_E->add(hit->E);
          mul++;
          //p = true
        }//addPadHit

        static bool GammaGate(double E, double Emin, double Emax){
          return Emin <= E && E <= Emax;
        }

        //code that runs after all events have been analyzed
      void terminate() override { //after analyzing we need to terminate the program
      AbstractSortedAnalyzer::terminate();
      gDirectory->WriteTObject(tree); //writes the output of the analysis to the tree in a file
      }//terminate

      void clear() {
        //We need to clear the memory - done manually in cpp
        //must clear all assigned variables used in analysis before going again
        mul = 0;
        p = g1 = g2 = pg1 = pg2 = b = bg = false;
        Eg1, Eg2, peakval, Egated= NAN;
        AUSA::clear(
        *v_id,
        *v_dir, *v_pos,
        *v_theta, *v_phi, *v_angle,
        *v_Edep, *v_fEdep, *v_bEdep,
        *v_FI, *v_BI, *v_FE, *v_BE, *v_FT, *v_BT,
        *v_E, *v_Ea
            );
        hits.clear();
        }//clear







TVector3 origin; // origin of the projectiles to be analyzed
//parameters from the constructor
shared_ptr<Setup> setupSpecs;
shared_ptr<Target> target;
double implantation_depth;
string isotopetype;
bool exclude_hpges, include_DSSSD_rim, include_spurious_zone, include_banana_cuts, include_beta_region;
Bool_t p, g1, g2, pg1, pg2, b, bg;
TelescopeTabulation *pU1P1, *pU2P2, *pU3P3, *pU4P4, *pU6P6;


unordered_set<Detector_frib *> detectors;
Detector_frib *U1, *U2, *U3, *U4, *U5, *U6, *P1, *P2, *P3, *P4, *P5, *P6, *G1, *G2;//, *U5, *U6, *P5, *P6

TTree *tree; //defines the tree in which we save the new parameters
int NUM;

UInt_t mul{}, CLOCK{}; //TPATTERN{}, TPROTONS{},
SortedSignal clock; //tpattern, tprotons, are these tprotons the time related to the measurements?

Double_t Eg1, Eg2, peakval, Egated;
unique_ptr<DynamicBranchVector<unsigned short>> v_id;
unique_ptr<DynamicBranchVector<TVector3>> v_dir, v_pos;
unique_ptr<DynamicBranchVector<double>> v_theta, v_phi, v_angle;
unique_ptr<DynamicBranchVector<double>> v_Edep, v_fEdep, v_bEdep;
unique_ptr<DynamicBranchVector<unsigned short>> v_FI, v_BI;
unique_ptr<DynamicBranchVector<double>> v_FE, v_BE;
unique_ptr<DynamicBranchVector<double>> v_FT, v_BT;
unique_ptr<DynamicBranchVector<double>> v_E, v_Ea;

vector<Hit> hits;
unique_ptr<EnergyLossRangeInverter> pSiCalc;
unique_ptr<EnergyLossRangeInverter> pAlCalc;
unique_ptr<EnergyLossRangeInverter> aSiCalc;
vector<unique_ptr<EnergyLossRangeInverter>> pTargetCalcs;
vector<unique_ptr<EnergyLossRangeInverter>> aTargetCalcs;

}; //Class GeneralAnalysis
#endif