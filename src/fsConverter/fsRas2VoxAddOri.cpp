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

string appName("this app converts surfs in RAS coordindate (lh.white ) to vox coordindate (paraview), and then adds the origin for paraview overlaid visualization");
string appExample("appName lh.white filled.mhd lh.white.vtk \n\
the lh.white is input, filled.mhd contains head information of volume, and lh.white.vtk is the output");
Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
Option<bool> vtkFormatOpt("-v",false,"if set, then input lh.white is in vtk format.", false,no_argument);


int main(int argc, char **argv) {
  OptionParser cmdParser(appName,appExample);
  cmdParser.add(helpOpt);
  cmdParser.add(vtkFormatOpt);

  if(4 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,3);

  CTriSurface mySurf(argv[1],!vtkFormatOpt.set());
  mySurf.Ras2VoxFS();
  mySurf.AddOriginFS(argv[2]);
  mySurf.SaveAs(argv[3]);
 
  return 0;
} //end of main function;

