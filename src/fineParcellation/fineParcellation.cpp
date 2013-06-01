/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012. 
 *	License : GPL v3 
 *  
 */

#include <iostream>
#include "kaimingCommon.h"
#include "options.h"
#include "fineParcellator.h"
#include "tictoc.hpp"

#define __KML__NEEDTIMECOST

using namespace Utilities;
using namespace KML;
using namespace std;

string appName("this app parcellates the cortex using the fMRI signals using Metis");
string appExample("appName surf.central.vtk fMRI.lowres.nii.gz fmri.t1.mhd output.vtk [options] \n\
the input surf.central.vtk is the central surface reconstructed using white and pial surfaces; \n\
fmri.lowres.nii.gz is the preprocessed fmri signals, registered to a lowres t1 image; \n\
the fmri.t1.mhd is the header information of lowres t1 image.\n\
minCut is a minCut without pca;\n\
maxCut is the maximium cut number; \n\
output is a vtk file with parcels.");

Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false);
Option<int> numCutOpt("-n",200, "the number of partitions, default: 200", false,requires_argument); 
Option<float> balanceOpt("-b",2000, "the balance of partitions, default: 2000", false, requires_argument); 
Option<float> minCorrelOpt("-c", 0.3, "the minimal correlation necessary to be considered as an edge, default: 0.3", false,requires_argument); 
Option<int> solutionNumberOpt("-s",50,"the number of metis solutions to choose from: default 50", false,requires_argument);
Option<bool> isInteractiveOpt("-i", false,"whether or not interactive, default: false", false); 
Option<int> simTypeOpt("-t",0,"the type of similarity, default: static correlation (0); valid values are 0, 1 ( minimun of sliding windows correlation)", false, requires_argument);
Option<int> windLengthOpt("-w",40,"the length of sliding window, only takes effect if simTypeOpt is set to 1, default: 40", false, requires_argument);
Option<string> regSphereOpt("-R","","the name of reg sphere, if set, the program will save partition information on the registered sphere.", false, requires_argument);
Option<int> maxBadMappingOpt("-M",10,"the max number of vertices that have bad mapping, i.e., can not find corresponding voxels, default: 10.", false, requires_argument);
Option<string> cortexLabelOpt("-l", "", "the cortex label file from freesurfer, only the cortex vertices will be used for graph partition. ", false, requires_argument); 
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
  cmdParser.add(regSphereOpt);
  cmdParser.add(maxBadMappingOpt);
  cmdParser.add(cortexLabelOpt); 
  cmdParser.add(checkCPVidOpt);
  cmdParser.add(usePCA4BOLDOpt); 
  cmdParser.add(outputPartsVtkOpt);
  cmdParser.add(randomParcelOpt);
  
 
  if(5 > argc ){
    cmdParser.usage(); 
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,4);


 

  CFineParcellator myParcellator ( argv[1],argv[2],argv[3]); 
  myParcellator.SetSimType(simTypeOpt.value());
  myParcellator.SetBadMappingThresh(maxBadMappingOpt.value());  
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
  if(cortexLabelOpt.set())
  {
          myParcellator.SetROICortexVertices(cortexLabelOpt.value());
  }    
   if(regSphereOpt.set())
  {
    myParcellator.SetRegName(regSphereOpt.value()); 
  }
  
  myParcellator.DoNcut4WholeBrain(numCutOpt.value()); 
  
  myParcellator.SaveNcutPartitions(outputPartsVtkOpt.set(),checkCPVidOpt.value());
    
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

