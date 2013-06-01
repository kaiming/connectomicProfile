/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012. 
 *	License : GPL v3 
 *  
 */

#include <iostream>
#include "kaimingCommon.h"
#include "options.h"
#include "triSurface.h"

using namespace Utilities;
using namespace KML;
using namespace std;

string appName("this app generates the central cortical surface using freesurfer white and pial surfaces");
string appExample("appName lh.white lh.pial lh.central.vtk \n\
lh.white lh.pial are lh pairs from freesurfer outputs, lh.central.vtk is the output;");

Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
Option<bool> vtkFormatOpt("-v",false,"if set, then inputs are in vtk format.", false,no_argument);



int main(int argc, char **argv) {
  OptionParser cmdParser(appName,appExample);
  cmdParser.add(helpOpt);
  cmdParser.add(vtkFormatOpt);
 
  if(4 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,3);

  CTriSurface mySurfWM(argv[1],!vtkFormatOpt.set()),mySurfPial(argv[2],!vtkFormatOpt.set());
  
  if(mySurfPial.GetNumOfPoints() != mySurfWM.GetNumOfPoints())
  {
    cout<<"unmatched white and pial surface, please check your inputs."<<endl;
    exit(1); 
  }

 
  for (int idx = 0 ; idx < mySurfWM.GetNumOfPoints(); ++idx) {
    const VectorType& coord1 = mySurfPial.GetPointCoords(idx);
    const VectorType& coord2 = mySurfWM.GetPointCoords(idx); 
    
    VectorType midCoord = coord2+coord1; 
    midCoord*=0.5;      
    mySurfWM.ChangePointCoords(idx,midCoord);
 
    
  } // end of for::idx
 
  
  mySurfWM.SaveAs(argv[3]);

 
  return 0;
} //end of main function;
