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

    ///get the cut information, max cut, normalized cut value; 
    float ComputeNCut(graph_t *graph, int *where, int npart);
    int ComputeGraphCut(graph_t *graph, int *where);

    int GetMaxNbrPartitions(void);
    void BuildPartitionNbr(void);

    void PrintPartition(const KML::CPartition& currentPartition);
    
    ///do pca for each partition for eigen value and vector; 
    void DoPCA4Partition( KML::CPartition& currentPartition );
    void DoPCA4AllPartitions();
  
    //find the center of partitions;
    void FindPartitionCenters();
    //calculate the area of a partition;
    void CalculateAreas4Partitions();
 
    // build the partitions after metis; 
    void BuildPartitions(vector< int >& mapNID2PID, vector< int >& metisLabel4NodeID, vector< CPartition >& allPartitions);

    //get the closest vertex in mni native for the surface, the deformation maps between fs_LR space and origi
    void GetIDMappingViaDeformationField(string deformationFieldFile, vector< int >& idMaps, int numPoints);

    //get id maps for two surfaces at same space ;
    void MapIDs4TwoSurfsAtSameSpace(CTriSurface& surf1, CTriSurface& surf2, vector<int>& idMap1to2, vector<int>& idMap2to1);
    
    //save the partition information:: status (ratio, area, etc); 
    void SavePartitionStatus();

    //save the partition meta information for nDICCCOLs ( all other information of the partition, including sphere coord); 
    void SavePartitionMetaAndVisuals();
public:

    CFineParcellatorHCP(string surfName, string fmriName, string fsLabels, float tr=0.72,bool isSurfAndVolSameConvention=false);
    ~CFineParcellatorHCP();
 
    //the file maps a caret sphere to fs std sphere;
    //regSphereName is the sphere.reg in caret space; deformationField is the projection between fs_LR_32K and MNI native; fsSphereName is the original fs surf at level 6;
    void ConvertCaretStd2FreeSurferStd(string regSphereName, string fsSphereName, string deformNameFsLR2Native, string deformNameNative2FsLR);
  
    //save vertex mapping between fs lv6 and mni std native; anatomy is "lh" or "rh"
    void SaveIDMapsHFsLv6Native(string anatomy="lh");
    
    
    //ncut the whole brain using kway partition;
    void PartitioningWholeBrainWithMetis(int nParts=200, int numSolutions=50, float balance=1000, int niter=10, float corrThresh=0.35 );
    
    //output all; 
    void SavePartitions(string outputBaseName); 
 
    //whether or not to use pca for partion bold signal;
    void SetUsePCA4Bold(bool usePCA);
 
    //if set (value == true); then will assign edges to random value within range [ 0, 100];
    void SetRandomParcellation(bool value);
   

} ;


}//end of namespace KML;

#endif
