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

   // This check is crucial and should be enabled
   if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
       return true; // Intersection is valid and within bounds
   }

   return false; // Intersection is outside the bounds
}

bool IsFrameHit(const TVector3& source, const TVector3& direction, 
                double ringRadius, double ringThickness) {
    
    // Ring is centered at (0,0,0) with radius ringRadius
    // Ring extends from Z = -ringThickness/2 to Z = +ringThickness/2
    
    double sx = source.X(), sy = source.Y(), sz = source.Z();
    double dx = direction.X(), dy = direction.Y(), dz = direction.Z();
    
    double halfThickness = ringThickness / 2.0;
    double ringBottomZ = -halfThickness;  // -0.3mm
    double ringTopZ = halfThickness;     // +0.3mm
    
    // We need to find where the ray intersects the cylindrical surface at radius = ringRadius
    // Ray equation: P(t) = source + t * direction
    // Cylindrical surface: x² + y² = ringRadius²
    // Substituting: (sx + t*dx)² + (sy + t*dy)² = ringRadius²
    
    // Expand to quadratic: at² + bt + c = 0
    double a = dx*dx + dy*dy;
    double b = 2.0 * (sx*dx + sy*dy);
    double c = sx*sx + sy*sy - ringRadius*ringRadius;
    
    // Special case: ray is parallel to Z-axis (dx = dy = 0)
    if (a < 1e-12) {
        // Ray has constant radial distance
        double constantRadius = sqrt(sx*sx + sy*sy);
        double tolerance = 0.1; // mm - adjust as needed for ring thickness
        
        if (abs(constantRadius - ringRadius) > tolerance) {
            return false; // Ray never comes close to ring radius
        }
        
        // Ray is at the right radius, check if it passes through Z-slab
        if (abs(dz) < 1e-12) {
            // Ray is completely parallel to XY plane
            return (sz >= ringBottomZ && sz <= ringTopZ);
        }
        
        // Ray is parallel to Z-axis, check if it intersects Z-slab
        double t_bottom = (ringBottomZ - sz) / dz;
        double t_top = (ringTopZ - sz) / dz;
        
        // Ray intersects Z-slab if either intersection is in the future (t >= 0)
        // or if source is already within the slab
        return (sz >= ringBottomZ && sz <= ringTopZ) || 
               (t_bottom >= 0) || (t_top >= 0);
    }
    
    // General case: solve quadratic equation
    double discriminant = b*b - 4.0*a*c;
    if (discriminant < 0) {
        return false; // No intersection with cylindrical surface
    }
    
    // Two intersection points with the cylindrical surface
    double sqrt_disc = sqrt(discriminant);
    double t1 = (-b - sqrt_disc) / (2.0*a);
    double t2 = (-b + sqrt_disc) / (2.0*a);
    
    // Check if either intersection point is within the ring's Z-slab
    for (double t : {t1, t2}) {
        // Calculate Z-coordinate at this intersection
        double z_intersection = sz + t * dz;
        
        // Check if this intersection is within the ring's Z-slab
        if (z_intersection >= ringBottomZ && z_intersection <= ringTopZ) {
            return true; // Ray hits the ring material
        }
    }
    
    // Additional check: does the ray pass through the Z-slab at a radius close to ringRadius?
    // This handles cases where the ray might graze the ring
    if (abs(dz) > 1e-12) {
        // Find where ray enters and exits the Z-slab
        double t_enter_slab = (ringBottomZ - sz) / dz;
        double t_exit_slab = (ringTopZ - sz) / dz;
        
        if (t_enter_slab > t_exit_slab) {
            std::swap(t_enter_slab, t_exit_slab);
        }
        
        // Check radial distance at slab entry and exit points
        double tolerance = 0.1; // mm
        
        for (double t : {t_enter_slab, t_exit_slab}) {
            double x = sx + t * dx;
            double y = sy + t * dy;
            double r = sqrt(x*x + y*y);
            
            if (abs(r - ringRadius) <= tolerance) {
                return true;
            }
        }
        
        // Check if minimum radial distance occurs within the Z-slab
        // The minimum occurs at t = -b/(2a) if this is within the slab
        double t_min_radius = -b / (2.0*a);
        double z_min = sz + t_min_radius * dz;
        
        if (z_min >= ringBottomZ && z_min <= ringTopZ) {
            double x_min = sx + t_min_radius * dx;
            double y_min = sy + t_min_radius * dy;
            double r_min = sqrt(x_min*x_min + y_min*y_min);
            
            if (abs(r_min - ringRadius) <= tolerance) {
                return true;
            }
        }
    }
    
    return false; // No intersection with ring material
}


int main(){
    string setup_dir = getProjectRoot() + "/setup/";
    shared_ptr<Target> target = make_shared<Target>(JSON::readTargetFromJSON(setup_dir+"target.json"));
    auto setup = JSON::readSetupFromJSON(setup_dir+"setupEE.json");


    double implantation_depth = 17/1e6; //implantation for Si
    vector<double> radius = {1e-6}; // mm; radius of the desired emission circle
    vector<string> sides = {"front", "back"};
    
    
    vector<string> detectors = {"U1", "U2", "U3", "U6"}; // for padvetoed calibration
    //vector<double> peakenergies = {385.72,904.02,1843.18,2076.74,2217.45}; //kev - for padvetoed calibration 
    //vector<string> detectors = {"U5"}; // for padvetoed calibration
    //vector<double> peakenergies = {904.02,1843.18,2076.74,2217.45,3337.75,4089.18,5402.61}; //kev - for padvetoed calibration 
    vector<double> peakenergies = {4089.18};
    
    int N = 100000;
    double rframe = 6; // mm
    double tframe = 1; // mm



    unique_ptr<EnergyLossRangeInverter> pSiCalc = defaultRangeInverter("p", "Silicon");
    unique_ptr<EnergyLossRangeInverter> pAlCalc = defaultRangeInverter("p","Aluminum");
    vector<unique_ptr<EnergyLossRangeInverter>> pTargetCalcs;
    for (auto &layer: target->getLayers()) {
        pTargetCalcs.push_back(defaultRangeInverter(Ion::predefined("p"), layer.getMaterial()));
        //eloss of protons in target material
      }

    for(const string& d : detectors){
    auto det = setup->getDSSD(d);
    
    /*
    Here one should introduce a dynamic read of that specific detector's different thicknesses for energy correction
    */
    
    auto dfdl = AUSA::getFrontDeadLayer(*setup->getDSSD(d));



    TVector3 origin = target->getCenter() + (target->getThickness()/2. - implantation_depth)*target->getNormal(); // origin of the point source
    //cout << "target:  "<< target->getThickness() << endl; 
    TVector3 normal = det->getNormal().Unit();
    //define detector bounds of the inner most strips:
    
    /*Need to define here a loop over strip boundaries --------*/
    for(const string& side : sides){
        cout << "# Detector=" << d <<  "\t" << side << " strips 2..15" << endl;
    for(int  s=1; s<15; ++s){ // if i = 1 return strip 2 if i==14 return strip 15
    TVector3 bound1, bound2, bound3, bound4;
    //TVector3 bound1 = det->getContinuousPixelPosition(1.5,1.5); // position between pixel 1 and 2 for both front strip and back strip, should return the corner pos of the  
    //TVector3 bound2 = det->getContinuousPixelPosition(15.5,1.5); // position between pixel 15 and 16 & 1 and 2 for front and back strip respectively
    //TVector3 bound3 = det->getContinuousPixelPosition(15.5,15.5);
    //TVector3 bound4 = det->getContinuousPixelPosition(1.5,15.5);

    /*
    Define borders for each strip in a loop over strips in a loop over sides.
    fx. if the side is "front" and strip is 7 then the bounds returns the 4 corners of the strip and the loop then checks whether or not the particle hits this strip

    Im a bit unsure here whether or not the limits should be 1.5 and 15.5 instead of 0.5 and 16.5. 
    Using 16.5 and 0.5 a given strip E is calculated by allowing events being recorded in 1-16 on the oppsosite side. 
    fx the front strip 6 energy is an average of all the energy deposited in backstrip 1-16
    */
    if(side=="front"){
    bound1 = det->getContinuousPixelPosition(s+0.5,0.5); 
    bound2 = det->getContinuousPixelPosition(s+1+0.5,0.5); 
    bound3 = det->getContinuousPixelPosition(s+1+0.5,16.5);
    bound4 = det->getContinuousPixelPosition(s+0.5,16.5);    
    }
    else if(side=="back"){
    bound1 = det->getContinuousPixelPosition(0.5,s+0.5); 
    bound2 = det->getContinuousPixelPosition(16.5,s+0.5); 
    bound3 = det->getContinuousPixelPosition(16.5,s+1+0.5);
    bound4 = det->getContinuousPixelPosition(0.5,s+1+0.5);    
    }
    else cerr << "Error in sides" << endl;
    

    
                                                     
    TVector3 EmissionPoint, direction, source, intersection;

    for(const double& e : peakenergies){
    double Etot = 0;
    int counter =0;
        //cout << "peakE" << e << endl;
    for(size_t j =0; j<radius.size(); j++){
    for(int i=0; i<N; i++){

    TVector3 source = PointGenerator(origin, radius[j]); // randomly generated point within a circle
    TVector3 direction = DirectionGenerator(); // randomly generated direction from the emission point

    /*
    Do we hit the target frame?
    Relevant for U5 & U6
    */
    if(!IsFrameHit(source,direction, rframe, tframe)){
        //cout << "frame is not hit! yay  " << endl;
    /*
    What is the intersection point between the randomly generated direction vector and the detector plane
    */
    if (FindIntersection(source, direction, bound1, bound2, bound4, normal, intersection)){
            //cout << "intersection.X  " << intersection.X() << "   intersection.Y" << intersection.Y() << "    intersection.Z" << intersection.Z() << endl;
            //cout << "E  " << e << endl;
            double E = e;
            /* Debugging statements
            double loss = 0;
            double loss2 = 0;
            double dist = 0;
            */
            double angle = (intersection-source).Angle(-det->getNormal()); // angle with respect to emission location
            //cout << "Angle "<< TMath::RadToDeg()*abs(angle) << endl;
            auto &to = intersection;
            for (auto &intersection: target->getIntersections(source, to)) {
            auto &calc = pTargetCalcs[intersection.index];
            //    loss += calc->getTotalEnergyLoss(E, intersection.transversed);
                E -= calc->getTotalEnergyLoss(E, intersection.transversed);
            //    dist += intersection.transversed;
            }//forloop for energy correction in target
            //cout << "Eloss from target   " << loss << endl;
            //cout << "target material traversed   " << dist*1e6 << "nm" << endl;
            //front dead layer
            double cor_dfdl = dfdl/abs(cos(angle));
            //loss2 = pSiCalc->getTotalEnergyLoss(E, cor_dfdl);
            E -= pSiCalc->getTotalEnergyLoss(E, cor_dfdl);
            //cout << "E  " << E << endl;
            //cout << "Eloss from frontdeadlayer  " << loss2 << endl;
            Etot+=E;
            if(E!=0) counter++;
          
    
    }//does the source hit the infitie plane extended by the detector surface?
}//is the frame hit?
    }//for i in N
Etot/=counter;/* Sum up E's and divide by N?*/

cout << Etot << " ";
}//for rad in rads

}//for e in energies
cout << "#Strip " << s+1 << endl;
}//for i in strips 
}//for side in sides
}
return 0;
}//main


