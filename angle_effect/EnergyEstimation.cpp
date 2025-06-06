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

bool IsFrameHit(const TVector3& dir, double ringRadius, double ringThickness) {
    double vz = dir.Z();
    if (vz == 0) return false; // Parallel to ring plane → never intersects z-slab

    // s at which vector enters and exits ring slab in z
    double s1 = (-ringThickness / 2.0) / vz;
    double s2 = (+ringThickness / 2.0) / vz;
    if (s1 > s2) std::swap(s1, s2); // Ensure s1 < s2

    // Only consider forward-going vectors
    if (s2 < 0) return false;

    // Clamp s1 to positive
    if (s1 < 0) s1 = 0;

    // Pick midpoint s in this slab
    double smid = (s1 + s2) / 2.0;

    // Compute radius at this s
    double x = smid * dir.X();
    double y = smid * dir.Y();
    double r = std::sqrt(x*x + y*y);

    // Check if within ring radius ± tolerance
    double tolerance = ringThickness/2; // Half of 0.6 mm width
    return (r >= ringRadius - tolerance && r <= ringRadius + tolerance);
}


int main(){
    string setup_dir = getProjectRoot() + "/setup/";
    shared_ptr<Target> target = make_shared<Target>(JSON::readTargetFromJSON(setup_dir+"target.json"));
    auto setup = JSON::readSetupFromJSON(setup_dir+"setup.json");
    string output_dir = getProjectRoot() + "/simdata";
    std::unique_ptr<TFile> out = make_unique<TFile>((output_dir + "/angleeffect.root").c_str(), "recreate");


    double implantation_depth = 17/1e6; //implantation for Si
    vector<double> peakenergies = {385.72,904.02,1843.18,2076.74,2217.45}; //kev; 
    vector<double> radius = {0.1}; // mm; radius of the desired emission circle
    vector<string> sides = {"front", "back"};
    
    double Øframe = 6;
    double tframe = 0.6;



    unique_ptr<EnergyLossRangeInverter> pSiCalc = defaultRangeInverter("p", "Silicon");
    unique_ptr<EnergyLossRangeInverter> pAlCalc = defaultRangeInverter("p","Aluminum");
    vector<unique_ptr<EnergyLossRangeInverter>> pTargetCalcs;
    for (auto &layer: target->getLayers()) {
        pTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
        //eloss of protons in target material
      }

    int N = 100000;
    auto det = setup->getDSSD("U6");
    string d = "U6";
    /*
    Here one should introduce a dynamic read of that specific detector's different thicknesses for energy correction
    */
    
    auto dfdl = AUSA::getFrontDeadLayer(*setup->getDSSD(d));



    TVector3 origin(0,0,-0.3);// = target->getCenter() + (target->getThickness()/2. - implantation_depth)*target->getNormal(); // origin of the point source
    TVector3 normal = det->getNormal().Unit();
    //define detector bounds of the inner most strips:
    

    /*Need to define here a loop over strip boundaries --------*/
    for(const string& side : sides){
        cout << "#" << side << endl;
    for(int  s=1; s<15; ++s){ // if i = 1 return strip 2 if i==14 return strip 15
    TVector3 bound1, bound2, bound3, bound4;
    //TVector3 bound1 = det->getContinuousPixelPosition(1.5,1.5); // position between pixel 1 and 2 for both front strip and back strip, should return the corner pos of the  
    //TVector3 bound2 = det->getContinuousPixelPosition(15.5,1.5); // position between pixel 15 and 16 & 1 and 2 for front and back strip respectively
    //TVector3 bound3 = det->getContinuousPixelPosition(15.5,15.5);
    //TVector3 bound4 = det->getContinuousPixelPosition(1.5,15.5);

    if(side=="front"){
    bound1 = det->getContinuousPixelPosition(s+0.5,1.5); 
    bound2 = det->getContinuousPixelPosition(s+1+0.5,1.5); 
    bound3 = det->getContinuousPixelPosition(s+1+0.5,15.5);
    bound4 = det->getContinuousPixelPosition(s+0.5,15.5);    
    }
    else if(side=="back"){
    bound1 = det->getContinuousPixelPosition(1.5,s+0.5); 
    bound2 = det->getContinuousPixelPosition(15.5,s+0.5); 
    bound3 = det->getContinuousPixelPosition(15.5,s+1+0.5);
    bound4 = det->getContinuousPixelPosition(1.5,s+1+0.5);    
    }
    else cerr << "Error in sides" << endl;
    

    
                                                     
    TVector3 EmissionPoint, direction, source, intersection;

    for(const double& e : peakenergies){
    double Etot = 0;
    int counter =0;
        //cout << "peakE" << e << endl;
    for(size_t j =0; j<radius.size(); j++){
    for(int i=0; i<N; i++){

    TVector3 EmissionPoint = PointGenerator(origin, radius[j]); // randomly generated point within a circle
    TVector3 source = EmissionPoint;
    TVector3 direction = DirectionGenerator(); // randomly generated direction from the emission point

    if(!IsFrameHit(direction, Øframe, tframe)){
        //cout << "frame is not hit! yay  " << endl;

    if (FindIntersection(source, direction, bound1, bound2, bound4, normal, intersection)){
        if(IsDetectorHit(intersection, bound1,bound2,bound3,bound4)){
            //SUCCES
            //real angle analysis
            //cout << "E  " << e << endl;
            double E = e;
            //cout << "source: " << "   x" << source.X() << "   y" <<  source.Y() << "   z" << source.Z() << endl;
            //cout << "intersection: " <<  intersection.X() << intersection.Y()<<intersection.Z() << endl;
            //cout << "origin: " <<  "   x" <<origin.X() << "   y" <<origin.Y()<<"   z" <<origin.Z() << endl;
            double angle = (intersection-source).Angle(-det->getNormal()); // angle with respect to emission location
            //cout << "RealAngle "<< TMath::RadToDeg()*abs(RealAngle) << "   FakeAngle " << TMath::RadToDeg()*abs(FakeAngle) << endl;
            auto &to = intersection;
            for (auto &intersection: target->getIntersections(source, to)) {
            auto &calc = pTargetCalcs[intersection.index];
             E -= calc->getTotalEnergyLoss(E, intersection.transversed);
            }//forloop for energy correction in target

            //front dead layer
            double cor_dfdl = dfdl/abs(cos(angle));
            E -= pSiCalc->getTotalEnergyLoss(E, cor_dfdl);
            //cout << "E  " << E << endl;
            Etot+=E;
            if(E!=0) counter++;
        }//is the active detector area hit?
    }//does the source hit the infitie plane extended by the detector surface?
}//Is the frame hit?
}//for i in N
Etot/=counter;/* Sum up E's and divide by N?*/

cout << Etot << " ";
}//for rad in rads

}//for e in energies
cout << "#Strip " << s+1 << endl;
}//for i in strips 
}//for side in sides
return 0;
}//main


