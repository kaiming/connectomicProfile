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

string appName("this app converts surfs in fs format to vtkformat ");
string appExample("appName lh.white lh.vtk ");
Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);

int main(int argc, char **argv) {
  OptionParser cmdParser(appName,appExample);
  cmdParser.add(helpOpt);

  if(3 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,2);

  CTriSurface mySurf(argv[1],true);
  mySurf.SaveAs(argv[2]);

  return 0;
} //end of main function;

