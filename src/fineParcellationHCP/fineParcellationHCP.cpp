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
string appExample("appName \n\
\t/Disk3T/100307/MNINonLinear/fsaverage_LR32k/100307.L.midthickness.32k_fs_LR.surf.gii \n\ \t/Disk3T/100307/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_s2.atlasroi.L.32k_fs_LR.func.gii  \n\ \t/Disk3T/100307/MNINonLinear/fsaverage_LR32k/100307.L.aparc.a2009s.32k_fs_LR.label.gii \n\ \t/Disk3T/100307/MNINonLinear/Native/100307.L.sphere.reg.native.surf.gii  \n\
\t/Disk1T/bin/pkg2/lh.sphere.ico6.reg \n\
\t/Disk3T/100307/MNINonLinear/fsaverage_LR32k/32k_fs_LR2native.L.deform_map \n\ \t/Disk3T/100307/MNINonLinear/fsaverage_LR32k/native232k_fs_LR.L.deform_map \n\
\tsub.100307.200 \n\
\t-a lh \n\
NOTE: argv[1]: midSurf; argv[2]: fMRI time series; argv[3]:cortical label from freesurfer; argv[4]: mninative.sphere.reg; \n\
argv[5]: lh.sphere.ico6.reg; argv[7]: deform_map(FsLR2MniNative); argv[8]: deform_map(MniNative2FsLR); argv[9]: ouput basename; \n\
[options] ");

Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false);
Option<int> numCutOpt("-n",200, "the number of partitions, default: 200", false,requires_argument); 
Option<float> balanceOpt("-b",2000, "the balance of partitions, default: 2000", false, requires_argument); 
Option<float> minCorrelOpt("-c", 0.35, "the minimal correlation necessary to be considered as an edge, default: 0.35", false,requires_argument); 
Option<int> solutionNumberOpt("-s",50,"the number of metis solutions to choose from: default 50", false,requires_argument);
Option<bool> usePCA4BOLDOpt("-p",false,"whether to use pca for partition bold, default false, i.e., use average", false); 
Option<bool> randomParcelOpt("-r",false," if set, will parcel the surface randomly, default: false", false); 
Option<string> anatomyOpt("-a","lh","the anatomy (lh or rh) of the surface, default: lh", false,requires_argument); 
Option<int> numInterationOpt("-i",10,"the number of interation for each solution. Default: 10", false,requires_argument); 

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
  cmdParser.add(anatomyOpt);
  cmdParser.add(usePCA4BOLDOpt); 
  cmdParser.add(randomParcelOpt);
  cmdParser.add(numInterationOpt);
  
 
  if(9 > argc ){
    cmdParser.usage(); 
    exit(EXIT_FAILURE);
  }
  cmdParser.parse_command_line(argc,argv,8);

  
  //setup metis params; 

 CFineParcellatorHCP myParcellator(argv[1],argv[2],argv[3]); 
 myParcellator.PartitioningWholeBrainWithMetis(numCutOpt.value(),solutionNumberOpt.value(),balanceOpt.value(),numInterationOpt.value(), minCorrelOpt.value());
 myParcellator.ConvertCaretStd2FreeSurferStd(argv[4],argv[5],argv[6],argv[7]);
 myParcellator.SaveIDMapsHFsLv6Native(anatomyOpt.value()); 
 myParcellator.SavePartitions(argv[8]);

  #ifdef __KML__NEEDTIMECOST  
  timer.toc();
  cout<<"CPU Time (ms) for fine parcellation : "<<timer.totalTimeMilliSec()<<endl;
  #endif 
 
  return 0;
 
} //end of main function;

