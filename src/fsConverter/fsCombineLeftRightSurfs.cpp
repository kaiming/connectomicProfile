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

string appName("this app combines the left and right surfs into a whole brain surf");
string appExample("appName lh rh whole.vtk [options] \n\
note: lh rh can also be any two cortical surfaces."); 
Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
Option<bool> vtkFormatOpt("-v",false,"if set, then input lh and rh are in vtk format.", false,no_argument);
Option<bool> keepPointLabel("-l",false, "if set, the point label will be kept. default: no", false,no_argument);
Option<string> lhLabelNameOpt("-L","PointLabel","the name of point label for lh, only work when -l is set. default: PointLabel", false,requires_argument);
Option<string> rhLabelNameOpt("-R","PointLabel","the name of point label for rh, only work when -l is set. default: PointLabel", false,requires_argument);

int main(int argc, char **argv) {
  OptionParser cmdParser(appName,appExample);
  cmdParser.add(helpOpt);
  cmdParser.add(vtkFormatOpt);
  cmdParser.add(keepPointLabel);
  cmdParser.add(lhLabelNameOpt);
  cmdParser.add(rhLabelNameOpt);

  if(4 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,3);
  
  CTriSurface mySurfLh(argv[1],!vtkFormatOpt.set());
  CTriSurface mySurfRh(argv[2],!vtkFormatOpt.set());
  if(keepPointLabel.set())
  {
    mySurfLh.ReadPointLabel(lhLabelNameOpt.value());
    mySurfRh.ReadPointLabel(rhLabelNameOpt.value());    
    mySurfLh.CombineSurfWithPointLable(mySurfRh);    
  }
  else
    mySurfLh.CombineSurf(mySurfRh);
  
  mySurfLh.SaveAs(argv[3]);
 
 
  return 0;
} //end of main function;

