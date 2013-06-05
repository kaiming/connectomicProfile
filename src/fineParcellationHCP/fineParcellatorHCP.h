/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012.
 *	License : GPL v3
 *
 */

#ifndef __FineParcellatorHCP__H
#define __FineParcellatorHCP__H

#include <iostream>
#include "indexer.h"
#include <vector>
#include "newimageall.h"
#include "kaimingCommon.h"
#include "kmlGPMetis.h"
#include "kmlPartitions.h"
#include "triSurface.h"

namespace KML
{
  
class CPrivacyFineParcellatorHCP; 
class CFineParcellatorHCP : public CTriSurface
{

private:
    CPrivacyFineParcellatorHCP* ndata;

    
    CFineParcellatorHCP& operator=(const CFineParcellatorHCP&);
    CFineParcellatorHCP(const CFineParcellatorHCP&);
    
    int GetMaxNbrPartitions(void);
    void BuildPartitionNbr(void);
  
    //this function corrects the topology issue of partitions;
    void CorrectParcelTopology();
    //check the topology of partitions; 
    void CheckPartitionTopology();

    void PrintPartition(const KML::CPartition& currentPartition);
    void DoPCA4Partition( KML::CPartition& currentPartition );
    void DoPCA4AllPartitions();
    
    void GenerateColors4Partitions();
    
    void OutputPartitionSizeHistogram(std::string fileName, int minSize, int maxSize);
    void PlotPrincipleVector(std::vector<float>& src);
    float GetSimilarityTwoTimeSeries(vector<float>& t1, vector<float>& t2);
    void CheckZeroPrincipalVectors4Parts(void);
     
    //find the center of partitions;
    void FindPartitionCenters();
    //calculate the area of a partition; 
    void CalculateAreas4Partitions();     
    //check scatter problem and fix; 
    void CheckScattersAndFix();
    void BuildPartitions(vector< int >& mapNID2PID, vector< int >& metisLabel4NodeID, vector< CPartition >& allPartitions);


public:

    CFineParcellatorHCP(string surfName, string fmriName, string fsLabels, float tr=0.72,bool isSurfAndVolSameConvention=false);
    ~CFineParcellatorHCP();
 
 
    //the base name of output files;
    int  SetOutputBaseName(string name); 
    
    //the file maps a caret sphere to fs std sphere; 
    void ConvertCaretStd2FreeSurferStd(string caretSphere, string fsSphere, string outputMap="");    

    //min threshold for correlation; 
    void SetMinCorrelThresh(float minThresh);
    float GetMinCorrelThresh(void)const; 
    
    //write out the partitions;
    void VisualizePartitionsOnSurface(std::string fileName);
   
    //ncut the whole brain using kway partition; 
    void PartitioningWholeBrainWithMetis(int nParts=200, int numSolutions=50, float balance=1000, int niter=10 );
 
    //save the partition information;
    //if savefile was set, then all information would be saved to the file;  
    void PrintAllPartitionInfo(bool isInteractive = false, string saveFile = std::string(""));
  
    //save the partition meta information for nDICCCOLs;
    void SavePartitionMetaInforSphere(string sphereReg, string fileName);
 
    //read the cortex on Surface; 
    void SetROICortexVertices(std::string fileName);
    
    //save  partitions, if checkpid set, then will calculate the cp for that pid; 
    void SavePartitions(bool toSaveVtk,  int pid = -1 );
    
    //whether or not to use pca for partion bold signal; 
    void SetUsePCA4Bold(bool usePCA); 
    
    float ComputeNCut(graph_t *graph, int *where, int npart);
    
    int ComputeGraphCut(graph_t *graph, int *where);
    
    //if set (value == true); then will assign edges to random value within range [ 0, 100]; 
    void SetRandomParcellation(bool value); 
    
} ;


}//end of namespace KML;

#endif
