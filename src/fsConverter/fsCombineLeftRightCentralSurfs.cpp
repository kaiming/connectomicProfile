/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012. 
 *	License : GPL v3 
 *  
 */

#include <iostream>
#include "kaimingCommon.h"
#include "options.h"
#include "newimageall.h"
#include "newmatio.h"
#include "fslio.h"
#include "triSurface.h" 


using namespace Utilities;
using namespace KML;
using namespace std;

string appName("this app combines the left,right and central surfs into a whole brain surf");
string appExample("appName lh.vtk rh.vtk cs.vtk whole.vtk [options] "); 
Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
Option<bool> addLabel4RenderingOpt("-l",false, "whether or not add label for rendering purpose, default: no", false,no_argument);

int main(int argc, char **argv) {
  OptionParser cmdParser(appName,appExample);
  cmdParser.add(helpOpt);
  cmdParser.add(addLabel4RenderingOpt);

  if(5 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,4);
  
  CTriSurface mySurfLh(argv[1]);
  CTriSurface mySurfRh(argv[2]);
  CTriSurface mySurfCs(argv[3]); 
  
  mySurfLh.CombineSurf(mySurfRh);
  mySurfLh.CombineSurf(mySurfCs);
  
  if(addLabel4RenderingOpt.set())
  {
    vector<int> & labels = mySurfLh.GetAllPointLabel();
    labels.resize(mySurfLh.GetNumOfPoints());
    int szCs = mySurfCs.GetNumOfPoints();
    int szRh = mySurfRh.GetNumOfPoints();
    int szLh = mySurfLh.GetNumOfPoints() - szCs - szRh;
    
    for ( int idx = 0; idx < szLh ; ++idx )
    {
      labels [ idx] =1;
    } // end for loop::idx
    for ( int idx = szLh; idx < szLh+szRh ; ++idx )
    {
      labels[idx] = 2; 
    } // end for loop::idx
    for ( int idx = szLh+szRh; idx < mySurfLh.GetNumOfPoints() ; ++idx )
    {
      labels[idx] = 3; 
    } // end for loop::idx
    
  }
  
  mySurfLh.NeedPointsLabel(true);
  mySurfLh.SaveAs(argv[4]);
  
  return 0;
} //end of main function;

