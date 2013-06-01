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

string appName("this app first removed the added ori, and then converts surfs in vox coordindate to freesurfer origianl RAS coordindate");
string appExample("appName lh.white filled.mhd lh.white.vtk \n\
note: the lh.white was added with ori previously, which was specified in filled.mhd \n\
lh.white is fs format. use mri_convert first if in vtk format");
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
  cmdParser.parse_command_line(argc,argv,4);

  CTriSurface mySurf(argv[1],!vtkFormatOpt.set());
  mySurf.RemoveOriginFS(argv[2]);
  mySurf.Vox2Ras();
  mySurf.SaveAs(argv[3]);
 
  return 0;
} //end of main function;

