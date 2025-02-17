#ifndef DETECTOR_FRIB_H
#define DETECTOR_FRIB_H

#include "gCut.h"

#include <string>
#include <utility>

#include <gsl/gsl_matrix.h>

#include <ausa/setup/Setup.h>
#include <ausa/setup/SquareDSSD.h>
#include <ausa/setup/PadDetector.h>

#include <telescope/TelescopeTabulation.h>

using namespace std;
using namespace Telescope;

namespace Detectors_ {
  enum DetectorType {
  NoType = 0,
  DSSSD,
  Pad,
  HPGe,
  };
  enum CalibrationType {
    NoCalibration = 0,
    Proton,
    Alpha,
    Gamma
    };

  class Detector_frib {
    public:
      Detector_frib()=default; // allows an instantiazition without specifying anything

      Detector_frib(unsigned short _id, string _name, DetectorType _type, CalibrationType _calibration, const shared_ptr<AUSA::Setup>& _setup)
      : id(_id), name(_name), type(_type), calibration(_calibration), setup(_setup) {
        switch (type) {
          case DSSSD:
            //if the detector is a DSSSD then it goes into the setup folder and finds a DSSSD detector with a
            //matching name and returns the front dead layer of this detector
            frontside_DeadLayer = AUSA::getFrontDeadLayer(*setup->getDSSD(this->name));
            // perhaps the naming of the original makes it more consistent with regards to the AUSA commands
            backside_DeadLayer = AUSA::getBackDeadLayer(*setup->getDSSD(this->name)); // same for backside
            Detector_thickness = AUSA::getThickness(*setup->getDSSD(this->name)); // same for thickness
            break;
          case Pad:
            frontside_DeadLayer = AUSA::getFrontDeadLayer(*setup->getSingleSided(this->name));
            backside_DeadLayer = AUSA::getBackDeadLayer(*setup->getSingleSided(this->name));
            Detector_thickness = AUSA::getThickness(*setup->getSingleSided(this->name));
            break;
          case HPGe:
            //nothing needs to be done
            break;
          case NoType:
            cerr << "Warning the detector type of " << this->name << "not recognized, set as NoType" << endl;
            break;
      }//switch
      }//standard detector instance

      Detector_frib(unsigned short _id, string _name, DetectorType _type,
               CalibrationType _calibration, const shared_ptr<AUSA::Setup>& _setup, double BetaCut)
      : Detector_frib(_id, _name, _type, _calibration, _setup) {
        setBetaCut(BetaCut);
      } //detector instance with betacutoff

   //define getters & setters

   unsigned short getId() const {return id;};
   string getName() const {return name;};
   DetectorType getType() const {return type;};
   CalibrationType getCalibration() const {return calibration;};
   Detector_frib* getPartner() const {return partner;};
   double getThickness() const {return Detector_thickness;};
   double getFrontDeadLayer() const {return frontside_DeadLayer;};
   double getBackDeadLayer() const {return backside_DeadLayer;};
   double getBetaCut() const {return BetaCut;};
   //get the banana cut for a specific detector
   gCut* getBananaCut() const {return banana;};
   //get telescope tabulations
   vector<TelescopeTabulation*> getTelescopeTabulations() const { return tTabulations; };

    //Setters
   void setId(unsigned short id) {this->id = id;};
   void setName(string name) {this->name = std::move(name);};
   void setType(DetectorType type) {this->type = type;};
   void setCalibration(CalibrationType calibration) {this->calibration = calibration;};
   void setBetaCut(double betaCut) {this->BetaCut = betaCut;};
   void setPartner(Detector_frib* partner) {this->partner = partner; this->withPartner=true;}
   void setThickness(double thickness) {this->Detector_thickness = thickness;}
   void setFrontDeadLayer(double frontdeadlayer) {this->frontside_DeadLayer = frontdeadlayer;}
   void setBackDeadLayer(double backdeadlayer) {this->backside_DeadLayer = backdeadlayer;}
   //set banana cut
   void setBananaCut(gCut* cut) {this->banana = cut;}
   //add telescope tabulations
   void addTelescopeTabulation(TelescopeTabulation* tT) { this->tTabulations.emplace_back(tT); }

   // If i have to make the addition of a solid angle load in then it can go here


   bool hasPartner() const {return withPartner;}


  private:
    unsigned short id{};
    string name;
    DetectorType type{};
    CalibrationType calibration{};
    Detector_frib* partner{};
    shared_ptr<AUSA::Setup> setup;
    double frontside_DeadLayer{}, Detector_thickness{}, backside_DeadLayer{}; // = 0.0 by default
    double BetaCut{};
    bool withPartner = false;
    gCut* banana{};
    vector<TelescopeTabulation*> tTabulations;

    }; //detector class

  void makePartners(Detector_frib* det1, Detector_frib* det2){
    det1->setPartner(det2);
    det2->setPartner(det1);
  }
  } //detectors namespace


#endif