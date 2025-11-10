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

    TVector3 P0 = bound1;
    double denom = normal.Dot(direction);
    if (std::abs(denom) < 1e-6) {
        return false;
    }

    double t = normal.Dot(P0 - source) / denom;
    if (t < 0) {
        return false;
    }
    intersection = source + t * direction;

    TVector3 bound3 = bound2 + (bound4 - bound1); 

    TVector3 edge1 = bound2 - bound1;
    TVector3 edge2 = bound4 - bound1;


    TVector3 toIntersection = intersection - bound1;
   double u = toIntersection.Dot(edge1) / edge1.Mag2();
   double v = toIntersection.Dot(edge2) / edge2.Mag2();

   if (u >= 0 && u <= 1 && v >= 0 && v <= 1) {
       return true;
   }

   return false; // Intersection is outside the bounds
}

bool IsFrameHit(const TVector3& source, const TVector3& direction, 
                double ringRadius, double ringThickness) {

    double sx = source.X(), sy = source.Y(), sz = source.Z();
    double dx = direction.X(), dy = direction.Y(), dz = direction.Z();
    
    double halfThickness = ringThickness / 2.0;
    double ringBottomZ = -halfThickness;  
    double ringTopZ = halfThickness;     
    
    double a = dx*dx + dy*dy;
    double b = 2.0 * (sx*dx + sy*dy);
    double c = sx*sx + sy*sy - ringRadius*ringRadius;
    
    if (a < 1e-12) {
        double constantRadius = sqrt(sx*sx + sy*sy);
        double tolerance = 0.1;
        if (abs(constantRadius - ringRadius) > tolerance) {
            return false; 
        }
        
        if (abs(dz) < 1e-12) {

            return (sz >= ringBottomZ && sz <= ringTopZ);
        }
        
        double t_bottom = (ringBottomZ - sz) / dz;
        double t_top = (ringTopZ - sz) / dz;
        
        return (sz >= ringBottomZ && sz <= ringTopZ) || 
               (t_bottom >= 0) || (t_top >= 0);
    }
    
 
    double discriminant = b*b - 4.0*a*c;
    if (discriminant < 0) {
        return false; 
    }

    double sqrt_disc = sqrt(discriminant);
    double t1 = (-b - sqrt_disc) / (2.0*a);
    double t2 = (-b + sqrt_disc) / (2.0*a);
    for (double t : {t1, t2}) {
        double z_intersection = sz + t * dz;
        if (z_intersection >= ringBottomZ && z_intersection <= ringTopZ) {
            return true;
        }
    }
    if (abs(dz) > 1e-12) {
        double t_enter_slab = (ringBottomZ - sz) / dz;
        double t_exit_slab = (ringTopZ - sz) / dz;
        
        if (t_enter_slab > t_exit_slab) {
            std::swap(t_enter_slab, t_exit_slab);
        }
        double tolerance = 0.1; 
        
        for (double t : {t_enter_slab, t_exit_slab}) {
            double x = sx + t * dx;
            double y = sy + t * dy;
            double r = sqrt(x*x + y*y);
            
            if (abs(r - ringRadius) <= tolerance) {
                return true;
            }
        }
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
    return false; //IsFrameHit returns false if the targetframe is not hit
}


int main(int argc, char* argv[]){
    string setup_path = getProjectRoot() + "/setup/setup.json";
    string target_path = getProjectRoot() + "/setup/target.json";
    string ion = "p";
    bool include_dssd_rim = false, no_target = false;
    double x = 0., y = 0., z = 0.;
    double radius = 1e-6;
    vector<string> sides = {"front", "back"};
    int N = 100000; // minimum to obtain no variances in energy - greatly reduces the runtime of the script...
    double rframe = 6; // mm - extension of target frame
    double tframe = 1; // mm - extension of target frame


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
        ("R, radius", "Radius of circular beam spot on target foil (mm)", cxxopts::value<string>()->default_value("Point source, radius=1e-6 mm"))
        ("n, iterations", "# of iterations for each simulation ", cxxopts::value<string>()->default_value("1e6"))
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
    z = result.count("z") ? result["z"].as<double>() : z;
    setup_path = result.count("setup") ? result["setup"].as<string>() : setup_path;
    auto setup = JSON::readSetupFromJSON(setup_path);
    target_path = result.count("target") ? result["target"].as<string>() : target_path;
    shared_ptr<Target> target;
    radius = result.count("radius") ? result["radius"].as<double>() : 1e-6;
    N = result.count("iterations") ? result["iterations"].as<int>() : 100000;
    double implantation_depth = 0.;
    TVector3 origin;
    if (!no_target) {// if a target exists we do...
        target = make_shared<Target>(JSON::readTargetFromJSON(target_path));
        implantation_depth = result.count("implantation") ?
                            stod(result["implantation"].as<string>()) : 1e6 * target->getThickness() / 2.; // nm
        implantation_depth /= 1e6; // mm
        origin = target->getCenter()
            + TVector3{x, y, -(target->getThickness()/2. - implantation_depth)};
    } else {
        target_path = "NOT USING TARGET";
        origin = TVector3{x, y, z};
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
        << "# Setup:               " << setup_path << endl
        << "# Target:              " << target_path << endl
        << "# Source point:        " << "(" << origin.X() << ", " << origin.Y() << ", "
                                        << origin.Z() << ")"  << endl
        << "# Radius of source:    " << radius << endl
        << "# Number of iteraitons:" << N << endl

        << "# Energies (keV):      ";
    for (const auto& E : result["energy"].as<vector<double>>()) {
        cout << E << " ";
    }
    cout << endl;


    
    for (const auto& d : result["detector"].as<vector<string>>()) {



    unique_ptr<EnergyLossRangeInverter> SiCalc = defaultRangeInverter(ion, "Silicon");
    unique_ptr<EnergyLossRangeInverter> AlCalc = defaultRangeInverter(ion,"Aluminum");
    vector<unique_ptr<EnergyLossRangeInverter>> TargetCalcs;
    if (target) {
    for (auto &layer: target->getLayers()) {
        TargetCalcs.push_back(defaultRangeInverter(Ion::predefined(ion), layer.getMaterial()));
        //eloss of given particle in target material
      }}

    auto det = setup->getDSSD(d);
    
    auto dfdl = AUSA::getFrontDeadLayer(*setup->getDSSD(d));
    //TVector3 origin = target->getCenter() + (target->getThickness()/2. - implantation_depth)*target->getNormal(); // origin of the point source
    //cout << "target:  "<< target->getThickness() << endl; 
    TVector3 normal = det->getNormal().Unit();
    
    for(const string& side : sides){
        cout << "# Detector=" << d <<  "\t" << side << " strips 2..15" << endl;
    for(int  s=1; s<15; ++s){ // if i = 1 return strip 2 if i==14 return strip 15 //////////
    TVector3 bound1, bound2, bound3, bound4;

    /*
    Define borders for each strip in a loop over strips in a loop over sides.
    fx. if the side is "front" and strip is 7 then the bounds returns the 4 corners of the strip and the loop then checks whether or not the particle hits this strip

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

    for (const auto& e : result["energy"].as<vector<double>>()) {
    double Etot = 0;
    int counter =0;
        //cout << "peakE" << e << endl;
    for(int i=0; i<N; i++){

    TVector3 source = PointGenerator(origin, radius); // randomly generated point within a circle
    TVector3 direction = DirectionGenerator(); // randomly generated direction from the emission point

    /*
    Do we hit the target frame?
    Relevant for U5 & U6
    */
    if(!IsFrameHit(source,direction, rframe, tframe) || !target){
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
            if (target) {
            auto &to = intersection;
            for (auto &intersection: target->getIntersections(source, to)) {
            auto &calc = TargetCalcs[intersection.index];
            //    loss += calc->getTotalEnergyLoss(E, intersection.transversed);
                E -= calc->getTotalEnergyLoss(E, intersection.transversed);
            //    dist += intersection.transversed;
            }}//forloop for energy correction in target
            //cout << "Eloss from target   " << loss << endl;
            //cout << "target material traversed   " << dist*1e6 << "nm" << endl;
            //front dead layer
            double cor_dfdl = dfdl/abs(cos(angle));
            //loss2 = pSiCalc->getTotalEnergyLoss(E, cor_dfdl);
            E -= SiCalc->getTotalEnergyLoss(E, cor_dfdl);
            //cout << "E  " << E << endl;
            //cout << "Eloss from frontdeadlayer  " << loss2 << endl;
            Etot+=E;
            if(E!=0) counter++;
            //if(E!=0) cout << TMath::RadToDeg()*abs(angle) << endl;

    
    }//does the source hit the plane extended by the detector surface?
}//is the frame hit?
    }//for i in N
Etot/=counter;

cout << Etot << " ";

}//for e in energies

cout << "#Strip " << s+1 << endl;

}//for i in strips 
}//for side in sides
}//for d in detectors
return 0;
}//main


