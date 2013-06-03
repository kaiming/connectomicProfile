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

namespace KML
{
class CPrivacyFineParcellatorHCP;
 
class CFineParcellatorHCP
{

private:
    CPrivacyFineParcellatorHCP* data;
    CFineParcellatorHCP& operator=(const CFineParcellatorHCP&);
    CFineParcellatorHCP(const CFineParcellatorHCP&);
    
    int GetMaxNbrParcels(void);
    int GetMaxNbrPartitions(void);
    void BuildPartitionNbr(void);

    void BuildInitParcels();
    void BuildPartitionsAfterNCut();
    void BuildParcelRelationship();

    //this function corrects the topology issue of partitions;
    void CorrectParcelTopology();
    //check the topology of partitions; 
    void CheckPartitionTopology();

    //the weight = floor( 10*correl(i,j) );
    void PrintPartition(const KML::CPartition& currentPartition);
    void DoPCA4Partition( KML::CPartition& currentPartition );
     
    
    void GenerateColors4Partitions();
    void OutputPartitionSizeHistogram(std::string fileName, int minSize, int maxSize);
    void RebuildPartitionsAndNbr();
    void PlotPrincipleVector(std::vector<float>& src);
    float GetSimilarityTwoTimeSeries(vector<float>& t1, vector<float>& t2);
    void CheckZeroPrincipalVectors4Parts(void);
    float GetDistanceHParcelAndPartition(int parcelIdx, int partitionIdx);
    double GetDistanceH2Parcels(int parcelIdx1, int parcelIdx2);
    double GetDistanceH2TimeSeries(vector<float>& s1, vector<float>& s2);
    void DetermineMeanSimilarity(int numSamples); 
    
    //find the center of partitions;
    void FindPartitionCenters();
    //calculate the area of a partition; 
    void CalculateAreas4Partitions(); 
    
    //check scatter problem and fix; 
    void CheckScattersAndFix();


public:

    CFineParcellatorHCP(std::string surfName, std::string fmriName, std::string fsLabels);
    ~CFineParcellatorHCP();


    void MapBolds2Surf(void);

    const CIndexer& GetIndexer(void);
    void SetIndexer(const CIndexer& othIndexer);

    //get parcel time-series;
    std::vector<float>& GetVertexTimeSeries(int parcelIdx);

    //get the cluster idx for a vertex;
    int GetVertexParcelIdx(int pointID);

    //get vertex time-series;
    std::vector<float>& GetVertexTimeSeries(int vertexID);

    void SetTR(float tr);
    float GetTR(void);

    void SetBadMappingThresh(int maxNum);
 
    int  SetBaseName(string name); 
    void SetCaretSurf2FSReg(string name);
    int GetBadMappingThresh(void);

    void SetMinCorrelThresh(float minThresh);
    float GetMinCorrelThresh(void)const; 
    
    void SetNumsMetisSolutions(int nSolutions); 

    void SetBalance(float minPercentageofVarance);
    float GetBalance(void);

    int GetNumVertices(void);
    int GetNumParcels(void);
    int GetNumFaces(void);

    int GetSimType();
    void SetSimType(int type);

    //only works when sim type=1 or sliding window similarity ;
    int GetWindLength();
    void SetWindLength(int windLength);

    //writeout the partitions;
    void SavePartitionsOnSurface(std::string fileName);

    //prepare the inputs for Ncuts, using partitions/parcels; ;
    void PrepareGraph4Metis(graph_t* graph, const KML::CPartition& currentPartition);
    //this will do the whole surf instead of parcels; 
    void PrepareGraph4Metis(graph_t* graph);
     
    //prepare the inputs for Ncuts;
//     void PrepareGraph4Ncut(const KML::CPartition& currentPartition,Eigen::MatrixXd& sims);
//     void PrepareGraph4Ncut(const KML::CPartition& currentPartition,Eigen::SparseMatrix<float>& sims);
    
    //ncut the whole brain using kway partition; 
    void PartitioningWholeBrainWithMetis( int nParts = 200 );

    //bi cut the brain;
    void BiNcut4WholeBrain(int minCuts, int maxCuts);

    //save the graph to graclus format;
    void SaveGraph(graph_t* graph, string fileName);

    //prune small partition;
    void PruneSmallPartitions(int minSize, float minCorrel=0.7);

    //save the partition information;
    //if savefile was set, then all information would be saved to the file; s
    void PrintAllPartitionInfo(bool isInteractive = false, string saveFile = std::string(""));



    //save the partition meta information for nDICCCOLs;
    void SavePartitionMetaInforSphere(string sphereReg, string fileName);

    //calculate the criteria of information metic;
    void CalculateInfoCritera(float& part1, float& aic, float& bic, float& icl);

    //calculate the Silhouette of the partitions; 
    double CalculateSilhouette();
	
    //determine the number of clusters;
    void DoPartitioning(int numParts);
    
    //read the cortex on Surface; 
    void SetROICortexVertices(std::string fileName);
    
    //save Ncut partitions, if checkpid set, then will calculate the cp for that pid; 
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
