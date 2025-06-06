#include <TFile.h>
#include "projectutil.h"
#include <iostream>
#include <TH1D.h>
#include <ausa/eloss/Default.h>
#include "TRandom.h"
#include <ausa/eloss/EnergyLossCalculator.h>
#include <utility>
#include <ausa/json/IO.h>
#include <telescope/DetectorTabulation.h>
#include <ausa/setup/Setup.h>
#include <ausa/setup/SquareDSSD.h>
#include <ausa/setup/PadDetector.h>



using namespace std;
using namespace AUSA::EnergyLoss;
using namespace AUSA;
using namespace EUtil;
using namespace Telescope;

TRandom randGen;
/*
Perhaps change the random generators to take and produce TVector3's instead...
*/
/* 
Random point generator within a circle of a given x,y and r
*/

TVector3 PointGenerator(TVector3 &center, double r){
    double theta = randGen.Uniform(0, 2 * M_PI); //random point from 0 to 2 pi
    double radius = r*sqrt(randGen.Uniform(0,1));
    double cartesian_x = center.X()+radius*cos(theta);
    double cartesian_y = center.Y()+radius*sin(theta);
    return TVector3(cartesian_x,cartesian_y,center.Z());
}//PointGenerator
/* 
Random direction generator from a given point x,y,z
*/
TVector3 DirectionGenerator(){
    double phi = randGen.Uniform(0, 2*M_PI); // azimuthal angle
    double theta = acos(randGen.Uniform(-1,1)); // polar angle
    double cartesian_x = sin(theta)*cos(phi);
    double cartesian_y = sin(theta)*sin(phi);
    double cartesian_z = cos(theta);
    TVector3 dir(cartesian_x,cartesian_y,cartesian_z);
    return dir;
}//DirectionGenerator

bool IsDetectorHit(const TVector3 &point, const TVector3 &bound1, const TVector3 &bound2, const TVector3 &bound3, const TVector3 &bound4){
    TVector3 v1 = bound1-point;
    TVector3 v2 = bound2-point;
    TVector3 v3 = bound3-point;
    TVector3 v4 = bound4-point;

    TVector3 n1 = v1.Cross(v2).Unit();
    TVector3 n2 = v2.Cross(v3).Unit();
    TVector3 n3 = v3.Cross(v4).Unit();
    TVector3 n4 = v4.Cross(v1).Unit();

    return (n1.Dot(n2) > 0) && (n2.Dot(n3) > 0) && (n3.Dot(n4) > 0);
}//IsDetectorHit

bool FindIntersection(const TVector3& source, const TVector3& direction, const TVector3& bound1, const TVector3& bound2, 
            const TVector3& bound4, const TVector3& normal, TVector3& intersection) {
// Compute a point on the plane (bound1 is used as P0)
    TVector3 P0 = bound1;

    // Calculate the dot product of the normal and direction
    double denom = normal.Dot(direction);
    if (std::abs(denom) < 1e-6) {
    // Line is parallel to the plane
        return false;
    }

    // Compute the parameter t for the line equation
    double t = normal.Dot(P0 - source) / denom;

    // Compute the intersection point
    intersection = source + t * direction;

    // Check if the intersection point lies within the bounds of the plane
    // (Assuming the plane is a rectangle defined by bound1, bound2, bound3, bound4)
    TVector3 bound3 = bound2 + (bound4 - bound1); // Compute bound3 from the given bounds

    // Vectors along the edges of the rectangle
    TVector3 edge1 = bound2 - bound1;
    TVector3 edge2 = bound4 - bound1;

    // Vector from bound1 to the intersection point
    TVector3 toIntersection = intersection - bound1;

    // Project the intersection point onto the plane's local coordinate system
    double u = toIntersection.Dot(edge1) / edge1.Mag2();
    double v = toIntersection.Dot(edge2) / edge2.Mag2();

    return true; // Intersection is valid and within bounds
}
int main(){
    string setup_dir = getProjectRoot() + "/setup/";
    shared_ptr<Target> target = make_shared<Target>(JSON::readTargetFromJSON(setup_dir+"target.json"));
    auto setup = JSON::readSetupFromJSON(setup_dir+"setup.json");
    string output_dir = getProjectRoot() + "/simdata";
    std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/angleeffect.root").c_str(), "recreate");




    unique_ptr<EnergyLossRangeInverter> pSiCalc = defaultRangeInverter("p", "Silicon");
    unique_ptr<EnergyLossRangeInverter> pAlCalc = defaultRangeInverter("p","Aluminum");
    vector<unique_ptr<EnergyLossRangeInverter>> pTargetCalcs;
    for (auto &layer: target->getLayers()) {
        pTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
        //eloss of protons in target material
      }

    //some values that could be read in through the argument instead for hardcoded if one wants to do multiple analyses
    double implantation_depth = 17/1e6; //implantation for Si
    double Epeak = 4089.18; //kev; most intense peak in Si
    vector<double> radius = {1,3,5}; // mm; radius of the desired emission circle
    int N = 1000000;
    auto det = setup->getDSSD("U2");
    string d = "U2";
    string p = "P2";
    /*
    Here one should introduce a dynamic read of that specific detector's different thicknesses for energy correction
    */
    
    auto dfdl = AUSA::getFrontDeadLayer(*setup->getDSSD(d));
    auto dbdl = AUSA::getBackDeadLayer(*setup->getDSSD(d));
    auto dbct = AUSA::getBackContactThickness(*setup->getDSSD(d));
    auto at = AUSA::getThickness(*setup->getDSSD(d)) - dfdl - dbdl - dbct;
    auto pfdl = AUSA::getFrontDeadLayer(*setup->getSingleSided(p));
    auto pfct = AUSA::getFrontContactThickness(*setup->getSingleSided(p));
    //cout << "active thickness :" << at << endl;


    TVector3 origin(0,-2,-0.3);// = target->getCenter() + (target->getThickness()/2. - implantation_depth)*target->getNormal(); // origin of the point source
    TVector3 normal = det->getNormal().Unit();
    //define detector bounds of the inner most strips:
    TVector3 bound1 = det->getContinuousPixelPosition(1.5,1.5); // position between pixel 1 and 2 for both front strip and back strip, should return the corner pos of the  
    TVector3 bound2 = det->getContinuousPixelPosition(15.5,1.5); // position between pixel 
    TVector3 bound3 = det->getContinuousPixelPosition(15.5,15.5);
    TVector3 bound4 = det->getContinuousPixelPosition(1.5,15.5);

    /*
    TVector3 boun = det->getContinuousPixelPosition(0, 0);
    TVector3 bouns = det->getContinuousPixelPosition(1, 1);
    TVector3 bounss = det->getContinuousPixelPosition(2, 2);
    TVector3 bound = det->getPixelPosition(1, 1);

    cout << "   bound: " << "x "<< bound.X() << "   y: " << bound.Y() << "   z: " << bound.Z() << endl;

    cout << "   boun: " << "x "<< boun.X() << "   y: " << boun.Y() << "   z: " << boun.Z() << endl;
    cout << "   bouns: " << "x "<< bouns.X() << "   y: " << bouns.Y() << "   z: " << bouns.Z() << endl;
    cout << "   bounss: " << "x "<< bounss.X() << "   y: " << bounss.Y() << "   z: " << bounss.Z() << endl;
*/
/*
    cout << "   bound1: " << "x "<< bound1.X() << "   y: " << bound1.Y() << "   z: " << bound1.Z() << endl;
    cout << "   bound2: " << "x "<< bound2.X() << "   y: " << bound2.Y() << "   z: " << bound2.Z() << endl;
    cout << "   bound3: " << "x "<< bound3.X() << "   y: " << bound3.Y() << "   z: " << bound3.Z() << endl;
    cout << "   bound4: " << "x "<< bound4.X() << "   y: " << bound4.Y() << "   z: " << bound4.Z() << endl;
    //cout << "Normal: " << "x" << normal.X() << "   y" << normal.Y() << "   z" << normal.Z() << endl;
  */  
    out->cd();
    map<string, TH1D*> hists;
    vector<string> rads = {"S", "M", "L"};
    for(const string& rad: rads){
        string namer = "RealE"+rad;
        hists[namer] = new TH1D(namer.c_str(),namer.c_str(), 120, 4030, 4150);
        string namef = "FakeE"+rad;
        hists[namef] = new TH1D(namef.c_str(),namef.c_str(), 120, 4030, 4150);
        string namerang = "RealAngles"+rad;
        hists[namerang] = new TH1D(namerang.c_str(),namerang.c_str(), 50, 0, 50);
        string namefang = "FakeAngles"+rad;
        hists[namefang] = new TH1D(namefang.c_str(),namefang.c_str(), 50, 0, 50);
        string namediff = "AnglesDiff"+rad;
        hists[namediff] = new TH1D(namediff.c_str(),namediff.c_str(),20,-20,20);   
    }
                                                                                        
    TVector3 EmissionPoint, direction, source, intersection;


    for(size_t j =0; j<rads.size(); j++){
    for(int i=0; i<N; i++){

    TVector3 EmissionPoint = PointGenerator(origin, radius[j]); // randomly generated point within a circle
    TVector3 source = EmissionPoint;
    TVector3 direction = DirectionGenerator(); // randomly generated direction from the emission point
    //cout << "EmissionPoint: " << "   x" << EmissionPoint.X() << "   y" <<  EmissionPoint.Y() << "   z" << EmissionPoint.Z() << endl;
    //cout << "direction: " << "   x" << direction.X() << "   y" <<  direction.Y() << "   z" << direction.Z() << endl;
        
    if (FindIntersection(source, direction, bound1, bound2, bound4, normal, intersection)){
        //cout << "yay: " <<  intersection.X() << intersection.Y()<<intersection.Z() << endl;
        if(IsDetectorHit(intersection, bound1,bound2,bound3,bound4)){
            //SUCCES
            //real angle analysis
            double E = Epeak;
            //cout << "source: " << "   x" << source.X() << "   y" <<  source.Y() << "   z" << source.Z() << endl;
            //cout << "intersection: " <<  intersection.X() << intersection.Y()<<intersection.Z() << endl;
            //cout << "origin: " <<  "   x" <<origin.X() << "   y" <<origin.Y()<<"   z" <<origin.Z() << endl;
            double RealAngle = (intersection-source).Angle(-det->getNormal()); // angle with respect to emission location
            double FakeAngle = (intersection-origin).Angle(-det->getNormal()); // angle with respect to center of target
            //cout << "RealAngle "<< TMath::RadToDeg()*abs(RealAngle) << "   FakeAngle " << TMath::RadToDeg()*abs(FakeAngle) << endl;
            auto &to = intersection;
            for (auto &intersection: target->getIntersections(source, to)) {
            auto &calc = pTargetCalcs[intersection.index];
             E -= calc->getTotalEnergyLoss(E, intersection.transversed);
            }//forloop for energy correction in target
            //cout << "target" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            //front dead layer
            double Real_dfdl = dfdl/abs(cos(RealAngle));
            double Fake_dfdl = dfdl/abs(cos(FakeAngle));
            E -= pSiCalc->getTotalEnergyLoss(E, Real_dfdl);
            //cout << "frontdeadlayer" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            //active layer -- assumed to be a big dead layer
            double Real_at = at/abs(cos(RealAngle));
            double Fake_at = at/abs(cos(FakeAngle));
            E -= pSiCalc->getTotalEnergyLoss(E, Real_at);
            //cout << "active layer" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            //back dead layer
            double Real_dbdl = dbdl/abs(cos(RealAngle));
            double Fake_dbdl = dbdl/abs(cos(FakeAngle));
            E -= pSiCalc->getTotalEnergyLoss(E, Real_dbdl);
            //cout << "backdeadlayer" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            //back contact thickness
            double Real_dbct = dbct/abs(cos(RealAngle));
            double Fake_dbct = dbct/abs(cos(FakeAngle));
            E -= pAlCalc->getTotalEnergyLoss(E, Real_dbct);
            //cout << "back contact thickness" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            //front contact thickness of pad
            double Real_pfct = pfct/abs(cos(RealAngle));
            double Fake_pfct = pfct/abs(cos(FakeAngle));
            E -= pAlCalc->getTotalEnergyLoss(E, Real_pfct);
            //cout << "front contact thickness" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            //front dead layer of pad
            double Real_pfdl = pfdl/abs(cos(RealAngle));
            double Fake_pfdl = pfdl/abs(cos(FakeAngle));
            E -= pSiCalc->getTotalEnergyLoss(E, Real_pfdl);
            //cout << "frontdeadlayer of pad" << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

            // now E should be the actual deposited energy in the pad only width of the peak is due to the stochastic nature of collisions
            //now we correct back with the wrong and right angles to see what we get. 
             /// in theory we dont need to correct the right angles back, but for consistency i find it a good idea to make sure there are no errors

            double realE = E;
            double fakeE = E;

            realE += pSiCalc->getTotalEnergyCorrection(realE, Real_pfdl);
            fakeE += pSiCalc->getTotalEnergyCorrection(fakeE, Fake_pfdl);

            realE += pAlCalc->getTotalEnergyCorrection(realE, Real_pfct);
            fakeE += pAlCalc->getTotalEnergyCorrection(fakeE, Fake_pfct);

            realE += pAlCalc->getTotalEnergyCorrection(realE, Real_dbct);
            fakeE += pAlCalc->getTotalEnergyCorrection(fakeE, Fake_dbct);

            realE += pSiCalc->getTotalEnergyCorrection(realE, Real_dbdl);
            fakeE += pSiCalc->getTotalEnergyCorrection(fakeE, Fake_dbdl);
            
            //i dont want to treat the active area as a dead layer as the measured quantity is not dependent on the angle difference
            fakeE += pSiCalc->getTotalEnergyCorrection(realE, Real_at);
            realE += pSiCalc->getTotalEnergyCorrection(realE, Real_at);

            realE += pSiCalc->getTotalEnergyCorrection(realE, Real_dfdl);
            fakeE += pSiCalc->getTotalEnergyCorrection(fakeE, Fake_dfdl);

            auto &from = intersection;
            for (auto &intersection: target->getIntersections(from, source)) {
            auto &calc = pTargetCalcs[intersection.index];
             realE += calc->getTotalEnergyCorrection(realE, intersection.transversed);
             fakeE += calc->getTotalEnergyCorrection(fakeE, intersection.transversed);
            }//forloop for energy correction in target
            


            //fill histograms
            string realkey = "RealAngles" + rads[j];
            string fakekey = "FakeAngles" + rads[j];
            string rekey = "RealE" + rads[j];
            string fekey = "FakeE" + rads[j];
            string diffkey = "AnglesDiff" + rads[j];
            hists[rekey]->Fill(realE);
            hists[fekey]->Fill(fakeE);
            hists[realkey]->Fill(TMath::RadToDeg()*RealAngle);
            hists[fakekey]->Fill(TMath::RadToDeg()*FakeAngle);
            hists[diffkey]->Fill(TMath::RadToDeg()*RealAngle-TMath::RadToDeg()*FakeAngle); // diff between real and fake angle
            
            //cout << "RealE: " << RealE << "  FakeE: " << FakeE << endl;

        }//is the active detector area hit?
    }//does the source hit the infitie plane extended by the detector surface?
};//for i in N
}//for rad in rads
    for (const auto& [key, hist] : hists) {
        hist->Write(); // Write each histogram to the output file
    }
    out->Close();
    return 0;
}//main

/*
define the angle to the center of the target and define the angle to the point of emission
    - based on this do an Eloss calc for both angles and write to two different histograms -> this should showcase the difference between real angle usage and assumed angle usage

*/

