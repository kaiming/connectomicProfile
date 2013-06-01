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

string appName("this app converts surfs in RAS coordindate (lh.white) to vox coordindate (paraview)");
string appExample("appName lh.white lh.white.vtk ");
Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
Option<bool> vtkFormatOpt("-v",false,"if set, then input lh.white is in vtk format.", false,no_argument);

int main(int argc, char **argv) {
  OptionParser cmdParser(appName,appExample);
  cmdParser.add(helpOpt);
  cmdParser.add(vtkFormatOpt);

  if(3 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,2);

  CTriSurface mySurf(argv[1],!vtkFormatOpt.set());
  mySurf.Ras2VoxFS();
  mySurf.SaveAs(argv[2]);

  return 0;
} //end of main function;

