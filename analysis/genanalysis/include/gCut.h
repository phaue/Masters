#ifndef GCUT_H
#define GCUT_H

#include <iostream>
#include <string>
#include <TGraph.h>
#include <TROOT.h>
#include <TFile.h>

#include "Cut.h"

using namespace std;

class gCut: public Cut{//inherit from the cut class
  public:
    gCut(const string &input_filename, const string &name_of_cut, cut_logic cutLogic=include_region) : Cut(cutLogic){
      TFile f(input_filename.c_str());
      // find the cut from the name and saves it to graph
      cut_graph = (TGraph*) gROOT -> FindObject(name_of_cut.c_str());
      }

      bool isInside(double x, double y) const override{
        return cut_graph -> IsInside(x, y);// overrides the bool from the cut class and now checks
        //whether this found graph object has x and y inside its specified region.
        }

      void print() {
        for(int i=0; i<cut_graph->GetN(); i++){
          cout << "(" << cut_graph->GetPointX(i) << ", " << cut_graph->GetPointY(i) << ")" << endl;
          }
      }
  private:
    TGraph *cut_graph;
    };
    #endif
