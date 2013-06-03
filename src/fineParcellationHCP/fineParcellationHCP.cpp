/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012. 
 *	License : GPL v3 
 *  
 */

#include <iostream>
#include "kaimingCommon.h"
#include "options.h"
#include "fineParcellatorHCP.h"
#include "tictoc.hpp"

#define __KML__NEEDTIMECOST

using namespace Utilities;
using namespace KML;
using namespace std;

string appName("this app parcellates the cortex using the fMRI signals using Metis");
string appExample("appName central.surf.gii/vtk bold4vertices.gii/arma aparc.a2009s output.vtk [options] \n\
the input central.surf.gii/vtk is the central surface reconstructed using white and pial surfaces; \n\
bold4vertices.gii/arma is the fmri signals, in gifti format or arma matrix \n\
aparc.a2009s is a file that shows the anatomical structures for surface, 0 is for unknown/medial wall.\n\
output is a vtk file with partitions.");

Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false);
Option<int> numCutOpt("-n",200, "the number of partitions, default: 200", false,requires_argument); 
Option<float> balanceOpt("-b",2000, "the balance of partitions, default: 2000", false, requires_argument); 
Option<float> minCorrelOpt("-c", 0.3, "the minimal correlation necessary to be considered as an edge, default: 0.3", false,requires_argument); 
Option<int> solutionNumberOpt("-s",50,"the number of metis solutions to choose from: default 50", false,requires_argument);
Option<bool> isInteractiveOpt("-i", false,"whether or not interactive, default: false", false); 
Option<int> simTypeOpt("-t",0,"the type of similarity, default: static correlation (0); valid values are 0, 1 ( minimun of sliding windows correlation)", false, requires_argument);
Option<int> windLengthOpt("-w",40,"the length of sliding window, only takes effect if simTypeOpt is set to 1, default: 40", false, requires_argument);
Option<string> caret2fsOpt("-R","","a file that maps fs_32K_LR to fs standard registered sphere.", false, requires_argument);
Option<int> checkCPVidOpt("-v",-1,"the vid to check its connectivity pattern, default -1, i.e., no check for any vid", false,requires_argument); 
Option<bool> usePCA4BOLDOpt("-p",false,"whether to use pca for partition bold, default false, i.e., use average", false); 
Option<bool> outputPartsVtkOpt("-o",false," if set, will output the colorized vtk surface, default: false",false); 
Option<bool> randomParcelOpt("-r",false," if set, will parcel the surface randomly, default: false", false); 

int main(int argc, char **argv) {
  
  //the following 4 lines are for computational time of fine parcellation; 
  //put it where you want to start timing; 
  #ifdef __KML__NEEDTIMECOST
  tictoc timer;
  timer.tic();
  #endif 
  
  OptionParser cmdParser(appName,appExample); 
  cmdParser.add(helpOpt); 
  cmdParser.add(numCutOpt);
  cmdParser.add(balanceOpt);
  cmdParser.add(minCorrelOpt);
  cmdParser.add(solutionNumberOpt);
  cmdParser.add(isInteractiveOpt);
  cmdParser.add(simTypeOpt);
  cmdParser.add(windLengthOpt);
  cmdParser.add(caret2fsOpt);
  cmdParser.add(checkCPVidOpt);
  cmdParser.add(usePCA4BOLDOpt); 
  cmdParser.add(outputPartsVtkOpt);
  cmdParser.add(randomParcelOpt);
  
 
  if(5 > argc ){
    cmdParser.usage(); 
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,4);


 

  CFineParcellatorHCP myParcellator ( argv[1],argv[2],argv[3]); 
  myParcellator.SetSimType(simTypeOpt.value());
  myParcellator.SetBaseName(argv[4]);
  myParcellator.SetBalance(balanceOpt.value());
  myParcellator.SetMinCorrelThresh(minCorrelOpt.value());
  myParcellator.SetUsePCA4Bold(usePCA4BOLDOpt.value());
  myParcellator.SetNumsMetisSolutions(solutionNumberOpt.value());
  myParcellator.SetRandomParcellation(randomParcelOpt.value()); 
  
  if(simTypeOpt.value()==1)
  {
    myParcellator.SetWindLength(windLengthOpt.value());
  }
 
  if(caret2fsOpt.set())
  {
    myParcellator.SetCaretSurf2FSReg(caret2fsOpt.value()); 
  }
  
  myParcellator.PartitioningWholeBrainWithMetis(numCutOpt.value()); 
  
  myParcellator.SavePartitions(outputPartsVtkOpt.set(),checkCPVidOpt.value());
    
  //the following 4 lines are for computational time of fine parcellation; 
  //put it where you want to end timing; 
  #ifdef __KML__NEEDTIMECOST  
  timer.toc();
  cout<<"CPU Time (ms) for fine parcellation : "<<timer.totalTimeMilliSec()<<endl;
  #endif 
     
  
//   fstream inStrm;
//   KML::OpenReadStrmBinary(inStrm,"lh.parts.200.CenterMetaInfo");
//   CPartition tmpPartition; 
//   vector<CPartition> allParts; 
//   while(inStrm>>tmpPartition)
//   {
//    allParts.push_back(tmpPartition);     
//   }
//   cout<<allParts.size()<<endl;

  return 0;
 
} //end of main function;

