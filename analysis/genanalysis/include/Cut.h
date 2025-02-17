#ifndef CUT_H
#define CUT_H

#include <iostream>
#include <string>
#include <TROOT.h>

using namespace std;

enum cut_logic { //creates an enumeration of two different cases either its in or its out
  include_region,
  exclude_region
  };

class Cut{
  public:
    explicit Cut(cut_logic cutLogic) : cutLogic(cutLogic){}; // ensures that cutLogic is set to either 0 or 1
    virtual bool isInside(double x, double y) const=0; //function to be overridden in future scripts
    //The isInside criteria is specified somewhere later, this is just a helper function so to say

    bool isSatisfied(double x, double y) const{
      return cutLogic ==include_region ? isInside(x, y): !isInside(x, y);
    //If cutLogic is equal to include_region then it checks if the point x,y is inside the region by calling
     // isInside, if however cutLogic is equal to exclude_region then it checks if the point x,y is outside the region
      //method to find whether or not a point x,y lies within a certain region or not
    }

    protected:
      cut_logic cutLogic; // makes sure that this variable cannot be changed outside the scope of this class
      };

      #endif