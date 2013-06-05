#include "fineParcellatorHCP.h"
#include <vector>
#include "ColorSchemeRGB.h"
#include "kmlPartitions.h"
#include "kmPCA.h"
#include "tictoc.hpp"

using std::vector;
using std::endl;
using std::cout;
using namespace KML;

#define __KML__NEEDTIMECOST

namespace KML {


class CPrivacyFineParcellatorHCP
{

    friend class CFineParcellatorHCP;

private:

    vector<CPartition> allPartitions;

    //partition neighbourhood;
    vector<set<int> > nbrPartitions;

    //keep records of all partition lables;
    vector<int> partionLabel4Vertex;

    //the areas of partitions;
    vector<float> allPartitionAreas;

    //whether pca or average for partition signal;
    bool usePCA4PartitionTimeSeries;

    //pca 1st eigen vector's variance fraction;
    vector<float> allPartition1stEVecFractions;

    CColorSchemeRGB myColorScheme;

    int numPartitions;
    int maxNbrs4Partitions;

    string baseName4Output;

    vector<VectorType> mapCaret2FsSphereCoods;
    vector<int> mapCaret2FsSphereVertexIDs;

    //map nodeID to vertex ID;
    vector<int> mapNID2PID;
    graph_t* graph;
    vector<int> metisLabel4NodeID;  //keep record of the partition result from metis;


public:
    CPrivacyFineParcellatorHCP():usePCA4PartitionTimeSeries(false),numPartitions(200),baseName4Output(""), graph(NULL) {};

};
CFineParcellatorHCP::CFineParcellatorHCP(string surfName, string fmriName, string fsLabels, float tr, bool isSurfAndVolSameConvention): CTriSurface(surfName)
{
    ndata = new CPrivacyFineParcellatorHCP;
    this->ReadPointLabelGifti(fsLabels);
    this->ReadTimeSeries4Vertices(fmriName,isSurfAndVolSameConvention);
    this->SetTR(tr);
    this->FilteringBoldsAfni();
    this->BuildNeighbors();
    ndata->graph = this->BuildGraphViaBold4Metis(this->GetAllPointLabel(),&(ndata->mapNID2PID),0);
}
CFineParcellatorHCP::~CFineParcellatorHCP()
{
    if(ndata->graph)
        FreeGraph(&(ndata->graph));
    delete ndata;
}
void CFineParcellatorHCP::PartitioningWholeBrainWithMetis(int nParts, int numSolutions, float balance, int niter)
{
    //check value ranges;
    KML::CheckValueRange<int>(nParts,2,this->GetNumOfPoints());
    KML::CheckValueRange<int>(niter,1,10000);
    KML::CheckValueRange<int>(numSolutions,1,2000);
    KML::CheckValueRange<float>(balance,1.f, 1e10);

    ndata->numPartitions = nParts;

    params_t* params = new params_t;
    KmlInitParams(params,nParts);
    params->balance=balance;
    params->ncuts=numSolutions;
    params->niter=niter;

    KmlGraphPartitionMetis(ndata->graph,params, ndata->metisLabel4NodeID, false);
    INFO_message("normalized cut value: %.3f",this->ComputeNCut ( ndata->graph, &(ndata->metisLabel4NodeID)[0],nParts ));
    INFO_message("metis cut edges: %d",this->ComputeGraphCut ( ndata->graph, &(ndata->metisLabel4NodeID)[0]));

    //note the metisLabel4NodeID is nid, instead of pid; build partiion with results;
    this->BuildPartitions(ndata->mapNID2PID, ndata->metisLabel4NodeID, ndata->allPartitions);

    //build partition neighbors;
    this->BuildPartitionNbr();

    //calculate pca;
    this->DoPCA4AllPartitions();

    //calculate area;
    this->CalculateAreas4Partitions();

    //get center vertex;
    this->FindPartitionCenters();
 


    delete params ;
};
void CFineParcellatorHCP::BuildPartitions(vector< int >& mapNID2PID, vector< int >& metisLabel4NodeID, vector< CPartition >& allPartitions)
{
    allPartitions.resize(ndata->numPartitions);
    ndata->partionLabel4Vertex.resize(this->GetNumOfPoints(),-1);


    for ( int idxNode = 0; idxNode < metisLabel4NodeID.size() ; ++idxNode )
    {
        int label = metisLabel4NodeID[idxNode];
        int pid = mapNID2PID[idxNode]; // nodeid = idxNode;
        ndata->partionLabel4Vertex[pid] = label; // -1 for medial wall;
        allPartitions[label].partitionLabel = label;
        allPartitions[label].parcelIds.push_back(pid);

    } // end for loop::idxNode

    //update partition size;
    for ( int idxPart = 0; idxPart < allPartitions.size() ; ++idxPart )
    {
        allPartitions[idxPart].numParcels = allPartitions[idxPart].parcelIds.size();
    } // end for loop::idxPart

    this->ColorizeSurfWLabels(ndata->partionLabel4Vertex);
    this->SaveAs("tst.vtk");
}
void CFineParcellatorHCP::BuildPartitionNbr(void )
{
    ndata->nbrPartitions.clear();
    ndata->nbrPartitions.resize(ndata->numPartitions);

    for ( int idxPart = 0; idxPart < ndata->allPartitions.size() ; ++idxPart )
    {
        vector<int>& mems = ndata->allPartitions[idxPart].parcelIds;
        for ( int idxV = 0; idxV < mems.size() ; ++idxV )
        {
            int pid= mems[idxV];
            int pidLabel = ndata->partionLabel4Vertex[pid];
            list<int> nbrs ;
            this->GetNeighboringPointsOfPoint(pid,nbrs);
            for(list< int >::iterator it=nbrs.begin(); it!=nbrs.end(); ++it)
            {
                int nbrLabel = ndata->partionLabel4Vertex[*it];
                if(-1!=nbrLabel && nbrLabel!= pid)
                    ndata->nbrPartitions[pidLabel].insert(nbrLabel);
            }
        } // end for loop::idxV
    } // end for loop::idxPart

    //get max nbr partitions;
    ndata->maxNbrs4Partitions = 0;
    for ( int idxP = 0; idxP < ndata->nbrPartitions.size() ; ++idxP )
    {
        ndata->maxNbrs4Partitions   = ndata->maxNbrs4Partitions < ndata->nbrPartitions[idxP].size() ? ndata->nbrPartitions[idxP].size() : ndata->maxNbrs4Partitions;
    } // end for loop::idxP
    INFO_message("max number of nbr patitions is: %d", ndata->maxNbrs4Partitions);
}
void CFineParcellatorHCP::CalculateAreas4Partitions()
{

    this->WarpPointLabel2CellLabel(ndata->partionLabel4Vertex);
    ndata->allPartitionAreas.resize(ndata->numPartitions+1,0);

    for ( int idxCell = 0; idxCell < this->GetNumOfCells() ; ++idxCell ) {
        int cellLabel = this->GetCellLabel ( idxCell );
        if ( -1!=cellLabel ) {
            ndata->allPartitionAreas[cellLabel]+= this->GetCellArea ( idxCell );
        }
    } // end for loop::idxCell

    for ( int idx = 0; idx < ndata->allPartitions.size() ; ++idx ) {
        ndata->allPartitions[idx].area = ndata->allPartitionAreas[idx];
    } // end for loop::idx

//     cout<<ndata->allPartitionAreas<<endl;
//     cout<<KML::StatsMean(ndata->allPartitionAreas)<<endl;


}
void CFineParcellatorHCP::DoPCA4Partition(CPartition& currentPartition)
{
    currentPartition.principalVector.clear();
    currentPartition.principalVector.resize ( this->GetNumOfTimePoints(),0 );

    //get all time series, and put it in arma;
    set<int> allVertices;
    for ( int idx = 0; idx < currentPartition.numParcels; ++idx ) {
        int pid = currentPartition.parcelIds[idx];
        allVertices.insert ( pid );
    } // end for loop::idx_t

    //check the numbe of time series;
    if(0==allVertices.size())
        ERROR_exit("there's no vertex for current partition");
    if(1==allVertices.size())
    {
        WARNING_message("there's only one vertex for current partition, will use it's time sereis as PCA.");
        currentPartition.ratio=1.;
        float* ptrCol=this->GetTimeSeries4Vertex(*(allVertices.begin()));
        memcpy((void*)(&(currentPartition.principalVector[0])),(void*)(ptrCol),sizeof(float)*this->GetNumOfTimePoints());
    }

    //do pca;
    CPCA pca;
    vector<vector<float> > rawdata, basis;
    vector<float> eigValue;
    {
        rawdata.resize(allVertices.size(), vector<float>(this->GetNumOfTimePoints(),0));
        int count=0;
        for( set<int>::iterator it=allVertices.begin(); it!= allVertices.end(); it++, count++)
        {
            int pid = *it;
            float* ptrTimeseriesPid = this->GetTimeSeries4Vertex(pid);
            memcpy((void*)(&(rawdata[count][0])),(void*)(ptrTimeseriesPid),sizeof(float)*this->GetNumOfTimePoints());
        }
    }

    pca.SetData ( rawdata );
    pca.DoPCA ( basis,eigValue );

    //get fration;
    {
        float sumEvalues = KML::StatsSum ( eigValue );

        if ( 0!=sumEvalues ) {
            currentPartition.ratio = eigValue[0]/ ( sumEvalues+FLT_MIN );
        } else {      //for sigular pca;
            WARNING_message("singular pca operation for partition : %d", currentPartition.partitionLabel);
            currentPartition.ratio=1;
        }
    }

    //get time series;
    {
        if ( ndata->usePCA4PartitionTimeSeries ) {

            for ( int idx = 0; idx < this->GetNumOfTimePoints() ; ++idx ) {
                currentPartition.principalVector[idx]= basis[idx][0];
            } // end for loop::idx

        } else {

            arma::fvec averageTimeSeries ( this->GetNumOfTimePoints());
            averageTimeSeries.fill ( 0 );

            for (set< int >::iterator it= allVertices.begin(); it!= allVertices.end(); it++)
            {
                float* ptrCol = this->GetTimeSeries4Vertex(*it);
                arma::fvec tmpVec(ptrCol,this->GetNumOfTimePoints(),false);
                averageTimeSeries+=tmpVec;
            }
            averageTimeSeries/= ( float ( allVertices.size() ) );
            for ( int idx = 0; idx < this->GetNumOfTimePoints() ; ++idx ) {
                currentPartition.principalVector[idx]=averageTimeSeries[idx];
            } // end for loop::idx
        }// done average;

    }; //end of get time series;

}
void CFineParcellatorHCP::DoPCA4AllPartitions()
{
  //start timing; 
  #ifdef __KML__NEEDTIMECOST
  tictoc timer;
  timer.tic();
  #endif 
  
      ndata->allPartition1stEVecFractions.resize(ndata->numPartitions,0);
    for ( int idx = 0; idx < ndata->numPartitions ; ++idx )
    {
        this->DoPCA4Partition(ndata->allPartitions[idx]);
        ndata->allPartition1stEVecFractions[idx]=ndata->allPartitions[idx].ratio;
    } // end for loop::idx
  
  //end timing; 
  #ifdef __KML__NEEDTIMECOST  
  timer.toc();
  INFO_message("CPU Time (ms) for DoPCA4AllPartitions : %d ms",timer.totalTimeMilliSec());
  #endif 
  

}
int CFineParcellatorHCP::ComputeGraphCut(graph_t* graph, int* where)
{
    int i, j, cut;

    if ( graph->adjwgt == NULL ) {
        for ( cut=0, i=0; i<graph->nvtxs; i++ ) {
            for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
                if ( where[i] != where[graph->adjncy[j]] )
                    cut++;
        }
    } else {
        for ( cut=0, i=0; i<graph->nvtxs; i++ ) {
            for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
                if ( where[i] != where[graph->adjncy[j]] )
                    cut += graph->adjwgt[j];
        }
    }

    return cut/2;

}
float CFineParcellatorHCP::ComputeNCut(graph_t* graph, int* where, int npart)
{
    int i, j, cm, nvtxs;
    idx_t *ncut, *degree, *xadj, *adjncy;
    float result;
    idx_t * adjwgt;

    ncut = ismalloc ( npart, 0, "ComputeNCut: ncut" );
    degree = ismalloc ( npart, 0, "ComputeNCut: degree" );
    nvtxs = graph->nvtxs;
    xadj = graph->xadj;
    adjncy = graph->adjncy;
    adjwgt = graph->adjwgt;

    if ( graph->adjwgt == NULL ) {
        for ( i=0; i<nvtxs; i++ ) {
            cm = where[i];
            for ( j=xadj[i]; j<xadj[i+1]; j++ ) {
                degree[cm] ++;
                if ( cm != where[adjncy[j]] )
                    ncut[cm] ++;
            }
        }
    } else {
        for ( i=0; i<nvtxs; i++ ) {
            cm = where[i];
            for ( j=xadj[i]; j<xadj[i+1]; j++ ) {
                degree[cm] += adjwgt[j];
                if ( cm != where[adjncy[j]] )
                    ncut[cm] += adjwgt[j];
            }
        }
    }
    int empty = 0;
    result =0;
    for ( i=0; i<npart; i++ ) {
        if ( degree[i] == 0 )
            empty++;
        if ( degree[i] >0 )
            result +=  ncut[i] *1.0/ degree[i];
    }
    //printf("Empty clusters: %d\n", empty);
    free ( ncut );
    free ( degree );
    return result+empty;
}
void CFineParcellatorHCP::FindPartitionCenters()
{
    // only do this after partiion is done;
    if(0==ndata->numPartitions)
        ERROR_exit("PartitioningWholeBrainWithMetis was not performed!");

    for ( int idxPart = 0; idxPart < ndata->allPartitions.size() ; ++idxPart ) {

        VectorType meanCoord;
        int count = ndata->allPartitions.operator[] ( idxPart ).parcelIds.size();
        for ( int idxPartVertices = 0; idxPartVertices < ndata->allPartitions.operator[] ( idxPart ).parcelIds.size() ; ++idxPartVertices ) {
            int vid = ndata->allPartitions.operator[] ( idxPart ).parcelIds[idxPartVertices];
            meanCoord += this->GetPointCoords ( vid );
        } // end for loop::idxPcls
        meanCoord/= count;

        //find the closest vertex;
        int closeVertexId = -1;
        float minDistance = 10000;

        for ( int idxPartVertices = 0; idxPartVertices < ndata->allPartitions.operator[] ( idxPart ).parcelIds.size() ; ++idxPartVertices ) {
            int vid = ndata->allPartitions.operator[] ( idxPart ).parcelIds[idxPartVertices];
            const VectorType& theCoord = this->GetPointCoords ( vid );
            VectorType disCoord = theCoord - meanCoord;
            if ( disCoord.NormSquare() < minDistance ) {
                closeVertexId  = vid;
                minDistance = disCoord.NormSquare();
            }
        } // end for loop::idxPcls

        if ( -1==closeVertexId ) {
            ERROR_exit("error while finding centerVertexID for partition: %d", ndata->allPartitions[idxPart].partitionLabel);
        } else {
            ndata->allPartitions[idxPart].centerVertexID = closeVertexId;
        }
    } // end for loop::idxPart
// cout for checking;
    for ( int idx = 0; idx < ndata->allPartitions.size() ; ++idx )
    {
        cout<<idx<<" "<<ndata->allPartitions[idx].centerVertexID<<endl;
    } // end for loop::idx

}






























};
