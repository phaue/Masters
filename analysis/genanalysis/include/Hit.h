#ifndef HIT_H
#define HIT_H

#include <TVector3.h>
#include <TLorentzVector.h>
#include "Detector_frib.h"

using namespace Detectors_;

//The structure of a hit is defined her with all the relevant proporties needed for analysis
struct Hit{
  Detector_frib* detector;

  unsigned short id; //what detector was hit?

  TVector3 direction, position, targetposition;
  double theta, phi, angle; //Angles associated with the hit, angle is the angle of incidence to the surface of the detector

  unsigned short FI, BI;
  double FE, BE;
  double FT, BT;
  double Edep; //average energy of FE and BE --> FE+BE /2

  double E, Ea, Eg; // energy? alpha energy and gamma energy.
  };

 #endif