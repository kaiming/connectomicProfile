/*
*	copyright by kaiming li (kaiming.li@emory.edu) 2012.
*	License : GPL v3
*
*/
#include <string>
#include <iostream>
#include "options.h"
#include "kaimingCommon.h" 
#include "triSurface.h"



using namespace Utilities;
using namespace KML;
using namespace std;
 

int main(int argc, char **argv) {

  string appName(argv[0]); 
  string appExample("removes the spikes vertices after affine transformation.\
  \n\nDetails:\n\
  \tthis program...\n\
  \te.g., appName inputSurf.vtk outputSurf.vtk \n"); 
  Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
  Option<float> spikeDistanceOpt("-d",10,"the max distance of an edge, default: 10mm",false,requires_argument);  
  
  
  OptionParser cmdParser(appName,appName+"\t"+appExample);
  cmdParser.add(helpOpt);
  cmdParser.add(spikeDistanceOpt);
 
  if(3 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
  }
 
  cmdParser.parse_command_line(argc,argv,2);
  CTriSurface mySurf(argv[1]);
  mySurf.BuildNeighbors();
//   mySurf.BuildEdge(); 
  
  set<int> allSpikeVertices; 
  for ( int idxFaces = 0; idxFaces < mySurf.GetNumOfCells() ; ++idxFaces )
  {
    const CellType& allVertices = mySurf.GetCell(idxFaces); 
    const VectorType& v1 = mySurf.GetPointCoords(allVertices.x);
    const VectorType& v2 = mySurf.GetPointCoords(allVertices.y);
    const VectorType& v3 = mySurf.GetPointCoords(allVertices.z);
    
    VectorType dis12 = v1-v2; 
    VectorType dis13 = v1-v3; 
    VectorType dis23 = v2-v3;

    int isLong12 = dis12.Norm() > spikeDistanceOpt.value() ? 1 : 0; 
    int isLong13 = dis13.Norm() > spikeDistanceOpt.value() ? 1 : 0; 
    int isLong23 = dis23.Norm() > spikeDistanceOpt.value() ? 1 : 0; 
    int totalValue = isLong12+ isLong13 + isLong23; 
    
    if(0==totalValue)
    {
      continue;
    }else if (2==totalValue) 
    {
      if(0==isLong12)
	allSpikeVertices.insert(allVertices.z); 
      else if ( 0 == isLong13)
	allSpikeVertices.insert(allVertices.y);
      else 
	allSpikeVertices.insert(allVertices.x); 
      
    }else if (3==totalValue || 1== totalValue ) 
    {
      cout<<"Error: spikeDistanceOpt is too small, I will exit."<<endl;
      exit(1);
    }else
    {
      cout<<"Unknow error: should not happen, will exit"<<endl;
      exit(2);      
    } 
  } // end for loop::idxFaces
 
  DEBUG_VAR(allSpikeVertices.size());
  if(allSpikeVertices.size())
  { 
    DEBUG_VAL("now removes spikes:"); 
    std::set< int >::iterator itAllSpikes = allSpikeVertices.begin(); 
    while(itAllSpikes != allSpikeVertices.end())
    {
      DEBUG_VAR(*itAllSpikes); 
      list<int> nbrs; 
      VectorType newCoords; 
      int currentSpikeVertex ; 
      float nbrCount(0); 
      {
	currentSpikeVertex = *itAllSpikes; 
	mySurf.GetNeighboringPointsOfPoint(currentSpikeVertex,nbrs);
	std::list< int >::iterator itNbrs = nbrs.begin(); 
	while(itNbrs != nbrs.end())
	{
	  if ( *itNbrs != currentSpikeVertex ) 
	  {
	    newCoords+= mySurf.GetPointCoords(*itNbrs); 	  
	    nbrCount+=1; 
	  }
	  itNbrs++; 
	}
	if(0!=nbrCount)
	{
	  newCoords/=nbrCount;
	  mySurf.SetPointComponent(currentSpikeVertex,0,newCoords.x);
	  mySurf.SetPointComponent(currentSpikeVertex,1,newCoords.y);
	  mySurf.SetPointComponent(currentSpikeVertex,2,newCoords.z);
	  
	}else
	{
	  cout<<"Error, no nbr points found! "<<__FUNCTION__<<" "<<__LINE__<<endl;; 
	  exit(1);
	}
	  
	
      }     
      
      itAllSpikes++;
    } // end for loop::idxSpikeVertex   
    
  }
 
 
  mySurf.SaveAs(argv[2]);
  
 
  return 0;
} //end of main function;


