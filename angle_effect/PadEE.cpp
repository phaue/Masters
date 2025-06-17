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
#include <cxxopts.hpp>




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

   // This check is crucial and should be enabled
   if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
       return true; // Intersection is valid and within bounds
   }

   return false; // Intersection is outside the bounds
}
vector<TVector3> findSquareCorners(const TVector3& center, double sideLength,
                                        const TVector3& normal, const TVector3& up) {
    double h = sideLength / 2.0;
    TVector3 v_up = up.Unit();
    TVector3 v_right = v_up.Cross(normal.Unit());

    TVector3 corner1 = center + (h * v_right) + (h * v_up);
    TVector3 corner2 = center - (h * v_right) + (h * v_up);
    TVector3 corner3 = center - (h * v_right) - (h * v_up);
    TVector3 corner4 = center + (h * v_right) - (h * v_up);

    return {corner1, corner2, corner3, corner4};
}


int main(){
    string setup_path = getProjectRoot() + "/setup/setup.json";
    string ion = "a";
    double x = 0., y = 0., z = 0.;
    double radius = 1e-6;
    int N = 100000; // minimum to obtain no variances in energy
    auto setup = JSON::readSetupFromJSON(setup_path);
    double detector_width = 50;
    vector<TVector3> source_origins = {
    {-5., 0., 0.},
    {-5., 0., 0.},
    {5., 0., 0.},
    {0., -5., 0.}
    };

    vector<string> detectors = {"P1", "P2", "P3", "P6"}; // for padvetoed calibration
    vector<double> peakenergies = {5155.4, 5485.74, 5804.96}; //kev - for padvetoed calibration 
    vector<TVector3> detector_centers = {
    {-31.9, -0.4, -31.9},
    {-31.9, -0.4, 31.9},
    {31.9, -0.4, 31.9},
    {0, -49, 0}
    };
    vector<TVector3> detector_normals = {
    {1, 0, 1},
    {1, 0, -1},
    {-1, 0, -1},
    {0, 1, 0}
    };
    vector<TVector3> detector_orientations = {
    {0, -1, 0},
    {0, -1, 0},
    {0, -1, 0},
    {1, 0, 1}
    };
    for (size_t i = 0; i< detectors.size(); ++i) {
        cout << "# " << detectors[i] << endl;
    unique_ptr<EnergyLossRangeInverter> SiCalc = defaultRangeInverter(ion, "Silicon");
    unique_ptr<EnergyLossRangeInverter> AlCalc = defaultRangeInverter(ion,"Aluminum");
    const TVector3& center = detector_centers[i];
    const TVector3& normal = detector_normals[i];
    const TVector3& orientation = detector_orientations[i];
    TVector3 origin = source_origins[i];
    double h = detector_width / 2.0;
    auto dfdl = AUSA::getFrontDeadLayer(*setup->getSingleSided(detectors[i]));
    auto dfct = AUSA::getFrontContactThickness(*setup->getSingleSided(detectors[i]));
    //TVector3 origin = target->getCenter() + (target->getThickness()/2. - implantation_depth)*target->getNormal(); // origin of the point source
    //cout << "target:  "<< target->getThickness() << endl; 
    TVector3 v_up = orientation.Unit();
    TVector3 v_right = v_up.Cross(normal.Unit());
    TVector3 bound1 = center + (h * v_right) + (h * v_up); // e.g., Top-Right
    TVector3 bound2 = center - (h * v_right) + (h * v_up); // e.g., Top-Left
    TVector3 bound3 = center - (h * v_right) - (h * v_up); // e.g., Bottom-Left
    TVector3 bound4 = center + (h * v_right) - (h * v_up); // e.g., Bottom-Right    

                                                 
    TVector3 EmissionPoint, direction, source, intersection;

    for (const auto& e : peakenergies) {
    double Etot = 0;
    int counter =0;
        //cout << "peakE" << e << endl;
    for(int i=0; i<N; i++){

    TVector3 source = PointGenerator(origin, radius); // randomly generated point within a circle
    TVector3 direction = DirectionGenerator(); // randomly generated direction from the emission point

    if (FindIntersection(source, direction, bound1, bound2, bound4, normal, intersection)){
            //cout << "intersection.X  " << intersection.X() << "   intersection.Y" << intersection.Y() << "    intersection.Z" << intersection.Z() << endl;
            //cout << "E  " << e << endl;
            double E = e;
            /* Debugging statements
            double loss = 0;
            double loss2 = 0;
            double dist = 0;
            */
            double angle = (intersection-source).Angle(-normal); // angle with respect to emission location
            double cor_dfct = dfct/abs(cos(angle));
            E -= SiCalc->getTotalEnergyLoss(E, cor_dfct);
            //front dead layer
            double cor_dfdl = dfdl/abs(cos(angle));
            //loss2 = pSiCalc->getTotalEnergyLoss(E, cor_dfdl);
            E -= SiCalc->getTotalEnergyLoss(E, cor_dfdl);
            //cout << "E  " << E << endl;
            //cout << "Eloss from frontdeadlayer  " << loss2 << endl;
            Etot+=E;
            if(E!=0) counter++;
            //if(E!=0) cout << TMath::RadToDeg()*abs(angle) << endl;

    
    }//does the source hit the infitie plane extended by the detector surface?
    }//for i in N
Etot/=counter;

cout << Etot << " ";

}//for e in energies
cout << endl;
}//for d in detectors
return 0;
}//main


