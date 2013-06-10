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


    //map between caret and fs; there's one map since they are at std space;FsLv6 :lh.sphere.ico6.reg, FsLR is the caret 32K surf;native is mni std native;  
    vector<int> mapIDNativeReg2FsLv6, mapIDFsLv6toNativeReg,mapIDFsLR2Native,mapIDNative2FsLR, mapIDFsLR2FsLv6,mapIDFsLv6toFsLR;

    // keep the sphere coordinate for all caret vertices;
    vector<VectorType> sphereCoordsAtFreeSurferStdSpace;

    bool doRandomParcellation;

public:
    CPrivacyFineParcellatorHCP():maxNbrs4Partitions(0),usePCA4PartitionTimeSeries(false),numPartitions(200),baseName4Output(""), graph(NULL),doRandomParcellation(false) {};


};
CFineParcellatorHCP::CFineParcellatorHCP(string surfName, string fmriName, string fsLabels, float tr, bool isSurfAndVolSameConvention): CTriSurface(surfName)
{
    ndata = new CPrivacyFineParcellatorHCP;
    this->ReadPointLabelGifti(fsLabels);
    this->ReadTimeSeries4Vertices(fmriName,isSurfAndVolSameConvention);
    this->SetTR(tr);
    this->FilteringBoldsAfni();
    this->BuildNeighbors();
}
CFineParcellatorHCP::~CFineParcellatorHCP()
{
    if(ndata->graph)
        FreeGraph(&(ndata->graph));
    delete ndata;
}
void CFineParcellatorHCP::PartitioningWholeBrainWithMetis(int nParts, int numSolutions, float balance, int niter, float corrThresh)
{
    //check value ranges;
    {
        KML::CheckValueRange<int>(nParts,2,this->GetNumOfPoints());
        KML::CheckValueRange<int>(niter,1,10000);
        KML::CheckValueRange<int>(numSolutions,1,2000);
        KML::CheckValueRange<float>(balance,1.f, 1e10);
	KML::CheckValueRange<float>(corrThresh,0.f,1);
    }
    ndata->numPartitions = nParts;
    params_t* params = new params_t;
    {
        KmlInitParams(params,nParts);
        params->balance=balance;
        params->ncuts=numSolutions;
        params->niter=niter;
    }

    ndata->graph = this->BuildGraphViaBold4Metis(this->GetAllPointLabel(),&(ndata->mapNID2PID),0,corrThresh);
    
    if(ndata->doRandomParcellation)
    {
        ERROR_exit("RandomParcellation not implimented for HCP data!");
    }

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

int CFineParcellatorHCP::GetMaxNbrPartitions(void )
{
    if(0==ndata->maxNbrs4Partitions)
        WARNING_message("partition nbrs has not been build.");

    return ndata->maxNbrs4Partitions;
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
void CFineParcellatorHCP::SetUsePCA4Bold(bool usePCA)
{
    ndata->usePCA4PartitionTimeSeries=usePCA;
}

void CFineParcellatorHCP::SetRandomParcellation(bool value)
{
    ndata->doRandomParcellation = value;
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
void CFineParcellatorHCP::PrintPartition(const CPartition& currentPartition)
{
    INFO_message("Partition info: \n\tlabel: %d\n\tnumofVertices: %d\n\tpca 1st eigen vector fraction: %.3f\n\tcenter vertex: %d\n\t%center vertex sphere coord: (%.3f,%.3f,%.3f)\n\tcortical area: %.3f\n\tnum of nbrs: %d", currentPartition.partitionLabel, currentPartition.numParcels,currentPartition.ratio,currentPartition.centerVertexID,currentPartition.sphereCoords.x,currentPartition.sphereCoords.y, currentPartition.sphereCoords.z,currentPartition.area,currentPartition.numNbrParts);

}
void CFineParcellatorHCP::SavePartitionStatus()
{
    vector<double> allRatios;
    vector<double> allSizes;
    arma::uvec histRatio,histSizes;
    {
        //set values for histgram;
        for ( int idx = 0; idx < ndata->numPartitions ; ++idx ) {
            allRatios.push_back ( ndata->allPartitions[idx].ratio );
            allSizes.push_back ( ndata->allPartitions[idx].parcelIds.size() );
        } // end for loop::idx
        arma::vec allRatiosArma ( &allRatios[0],allRatios.size(),false );
        arma::vec allSizesArma ( &allSizes[0],allSizes.size(),false );
        histRatio = arma::hist ( allRatiosArma,arma::linspace<arma::vec> ( 0,1,101 ) );
        histSizes = arma::hist ( allSizesArma,arma::linspace<arma::vec> ( 0,2000,101 ) );
    }

    {
        string fileName = ndata->baseName4Output+".partInfo"; 
        fstream outStrm;
        KML::OpenWriteStrmAscii ( outStrm,fileName );
        outStrm<< "number of partitions:\t"<< ndata->allPartitions.size() <<endl;
        outStrm<< "histogram of PCA ratio:"<<endl;
        outStrm<< histRatio.t() <<endl;
        outStrm<< "histogram of partition size(number of vertices):"<<endl;
        outStrm<< histSizes.t() <<endl;
        outStrm<< "areas of partitions: "<<endl;
        outStrm<< ndata->allPartitionAreas<<endl;
        outStrm.close();
    }

}
void CFineParcellatorHCP::SavePartitionMetaAndVisuals()
{
    string& baseName= ndata->baseName4Output;
    //get the r theta, phi; information;
    vector<VectorType> coordsSurfs;
    {
        for ( int idxPart = 0; idxPart < ndata->allPartitions.size() ; ++idxPart ) {
            const VectorType& coord = this->GetPointCoords ( ndata->allPartitions[idxPart].centerVertexID );
            coordsSurfs.push_back ( coord );
        } // end for loop::idxPart

        //write center as vertices; on native surface;
        string centerOnSurfName = baseName + ".CentersonSurf.vtk";
        KML::VisualizePointsUsingVerticeByVTK ( coordsSurfs, centerOnSurfName );
    }

    {
        //write center as vertices; on fs Lv6 sphere;
        string centerOnSphereName = baseName + ".CentersonSphere.vtk";
        KML::VisualizePointsUsingVerticeByVTK ( ndata->sphereCoordsAtFreeSurferStdSpace,centerOnSphereName );
    }
    //write metainfo
    fstream outStrm;
    string metaName = baseName+".CenterMetaInfo";
    {
        KML::OpenWriteStrmBinary ( outStrm,metaName );
        for ( int idxPart = 0; idxPart < ndata->allPartitions.size() ; ++idxPart ) {
            //write nbr partitions;
            ndata->allPartitions[idxPart].nbrPartitions.clear();
            ndata->allPartitions[idxPart].numNbrParts=ndata->nbrPartitions[idxPart].size();

            set<int>::iterator it = ndata->nbrPartitions[idxPart].begin();
            while ( it != ndata->nbrPartitions[idxPart].end() ) {
                ndata->allPartitions[idxPart].nbrPartitions.push_back ( *it );
                it++;
            }
            ndata->allPartitions[idxPart]>>outStrm;
        } // end for loop::idxPart
    }
    outStrm.close();
}

void CFineParcellatorHCP::SavePartitions(string outputBaseName)
{
  ///save partitions on surfs; 
  ndata->baseName4Output = outputBaseName; 
  {
    string fileName = ndata->baseName4Output+".vtk"; 
    this->ColorizeSurfWLabels(ndata->partionLabel4Vertex);
    this->SaveAs(fileName);
  }
  
  this->SavePartitionStatus();
  this->SavePartitionMetaAndVisuals();
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
//     for ( int idx = 0; idx < ndata->allPartitions.size() ; ++idx )
//     {
//         cout<<idx<<" "<<ndata->allPartitions[idx].centerVertexID<<endl;
//     } // end for loop::idx

}
void CFineParcellatorHCP::ConvertCaretStd2FreeSurferStd(string regSphereName, string fsSphereName, string deformNameFsLR2Native, string deformNameNative2FsLR )
{
    //start timing;
#ifdef __KML__NEEDTIMECOST
    tictoc timer;
    timer.tic();
#endif

    CTriSurface regSurf(regSphereName), fsLevel6Surf(fsSphereName);

    this->GetIDMappingViaDeformationField(deformNameFsLR2Native,ndata->mapIDNative2FsLR,regSurf.GetNumOfPoints()); // now we have the map fs_LR to mni native;
    this->GetIDMappingViaDeformationField(deformNameNative2FsLR,ndata->mapIDFsLR2Native, this->GetNumOfPoints()); // now we have the map fs_LR to mni native;
    this->MapIDs4TwoSurfsAtSameSpace(regSurf,fsLevel6Surf,ndata->mapIDNativeReg2FsLv6,ndata->mapIDFsLv6toNativeReg);

    //finally, we can map the fs_LR_id to fs_level6 sphere;
    ndata->mapIDFsLR2FsLv6.resize(this->GetNumOfPoints());
    ndata->mapIDFsLv6toFsLR.resize(fsLevel6Surf.GetNumOfPoints()); 
    ndata->sphereCoordsAtFreeSurferStdSpace.resize(this->GetNumOfPoints());
    for ( int idx = 0; idx < this->GetNumOfPoints() ; ++idx )
    {
        int nativeId= ndata->mapIDFsLR2Native[idx];
        int fsLv6Id = ndata->mapIDNativeReg2FsLv6[nativeId];
        ndata->mapIDFsLR2FsLv6[idx] = fsLv6Id;
        KML::ConvertCoordsEuclidean2Spherical(fsLevel6Surf.GetPointCoords(fsLv6Id),ndata->sphereCoordsAtFreeSurferStdSpace[idx]);
    } // end for loop::idx
    for ( int idx = 0; idx < fsLevel6Surf.GetNumOfPoints() ; ++idx )
    {
      int nativeId = ndata->mapIDFsLv6toNativeReg[idx];
      int fsLRId = ndata->mapIDNative2FsLR[nativeId]; 
      ndata->mapIDFsLv6toFsLR[idx] = fsLRId;       
    } // end for loop::idx

    //assign coordsSpheres to partitions;
    if(0==ndata->numPartitions)
        ERROR_exit("cortical parcellation not performed, please call PartitioningWholeBrainWithMetis() first!");

    for ( int idxP = 0; idxP < ndata->numPartitions ; ++idxP )
    {
        ndata->allPartitions[idxP].sphereCoords = ndata->sphereCoordsAtFreeSurferStdSpace[ndata->allPartitions[idxP].centerVertexID];
    } // end for loop::idxP

    //end timing;
#ifdef __KML__NEEDTIMECOST
    timer.toc();
    INFO_message("CPU Time for ConvertCaretStd2FreeSurferStd : %d ms",timer.totalTimeMilliSec());
#endif

}
void CFineParcellatorHCP::GetIDMappingViaDeformationField(string deformationFieldFile, vector< int >& idMaps, int numPoints)
{
    fstream deformStrm;
    KML::OpenReadStrmAscii(deformStrm,deformationFieldFile);

    string dataStartKEY("DATA-START");
    int numOfVertices=0;
    {   //get to the right location;
        string tmpString;
        while(getline(deformStrm,tmpString) && tmpString !=dataStartKEY)
        {};
        deformStrm>>numOfVertices;
        if(numOfVertices != numPoints)
            ERROR_exit("the number of vertices (%d) in deformation field (%s) is not equal to assumed (%d)", numOfVertices,deformationFieldFile.c_str(),numPoints);
    }

    //read data and determine the closing vertex;
    idMaps.resize(numOfVertices);
    for ( int idxV = 0; idxV < numOfVertices ; ++idxV )
    {
        int vid, node1,node2,node3,theNode(-1);
        float area1,area2,area3;
        deformStrm>>vid>>node1>>node2>>node3>>area1>>area2>>area3;
//         cout<<vid<<" "<<node1<<" "<<node2<<" "<<node3<<" "<<area1<<" "<<area2<<" "<<area3<<endl;

        //now determine which one is closest;
        VectorType p1(1,0,0),p2(0,1,0),p3(0,0,1);
        float allArea=area3+area2+area1;
        VectorType p = p1*(area1/allArea)+p2*(area2/allArea)+p3*(area3/allArea);
        VectorType disPP1(p-p1), disPP2(p-p2),disPP3(p-p3);
        float minDis = std::min(std::min(disPP1.Norm(),disPP2.Norm()),disPP3.Norm());
        if(minDis == disPP1.Norm())
            theNode = node1;
        else if (minDis == disPP2.Norm())
            theNode = node2;
        else
            theNode = node3;

        if(-1==theNode)
            ERROR_exit("can not find corresponding node or vertex: ",idxV);
        idMaps[idxV]=theNode;
    } // end for loop::idxV

}

void CFineParcellatorHCP::MapIDs4TwoSurfsAtSameSpace(CTriSurface& surf1, CTriSurface& surf2, vector< int >& mapID1to2, vector< int >& mapID2to1)
{

    CIndexer gridIndex;
    gridIndex.SetOffset(VectorType(-120,-120,-120));
    gridIndex.SetDims(VectorType(10,10,10));
    gridIndex.SetSize(Vector3D<size_t>(25,25,25));
    map<Vector3D<size_t>,vector<int> > idxCaret, idxFS;

    //build the mapping;
    {
        for ( int idx = 0; idx < surf1.GetNumOfPoints() ; ++idx )
        {
            const VectorType& coord = surf1.GetPointCoords(idx);
            Vector3D<size_t> grid = gridIndex.GetIndex(coord);
            idxCaret[grid].push_back(idx);
        } // end for loop::idx

        for ( int idx = 0; idx < surf2.GetNumOfPoints() ; ++idx )
        {
            const VectorType& coord = surf2.GetPointCoords(idx);
            Vector3D<size_t> grid = gridIndex.GetIndex(coord);
            idxFS[grid].push_back(idx);
        } // end for loop::idx
    }

    mapID1to2.resize(surf1.GetNumOfPoints());
    mapID2to1.resize(surf2.GetNumOfPoints());
    //coarse resolution mapping; might not be optimal; just good enough;

    for(map< Vector3D< size_t >, vector< int > >::iterator it=idxCaret.begin(); it!=idxCaret.end(); it++)
    {
        const Vector3D<size_t>& key = it->first;
        vector<int>& elemsCaret = it->second;
        vector<int>& elemsFS = idxFS[key];

        if(0==elemsFS.size()) // counterpart also has elems;
        {
            //search neighbors??? do whole brain , since there are very few;
            WARNING_message("cornered (%d) vertices, will search for whole brain ",elemsCaret.size());
            elemsFS.resize(surf2.GetNumOfPoints());
            for ( int idx = 0; idx < surf2.GetNumOfPoints() ; ++idx )
            {
                elemsFS[idx] = idx;
            } // end for loop::idx

        }

        for ( int idxElemCaret = 0; idxElemCaret < elemsCaret.size() ; ++idxElemCaret )
        {
            arma::fvec allDistance = arma::randu<arma::fvec>(elemsFS.size());
            const VectorType& coordCaret = surf1.GetPointCoords(elemsCaret[idxElemCaret]);
            for ( int idxElemFS = 0; idxElemFS < elemsFS.size() ; ++idxElemFS )
            {
                const VectorType& coordFS = surf2.GetPointCoords(elemsFS[idxElemFS]);
                VectorType dis = coordFS-coordCaret;
                allDistance[idxElemFS] = dis.Norm();
            } // end for loop::idxElemFS
            arma::uword minPos=0;
            float minDis = allDistance.min(minPos);
            mapID1to2[elemsCaret[idxElemCaret]] = elemsFS[minPos];
// 	    mapIDFS2Caret[elemsFS[minPos]]=elemsCaret[idxElemCaret];
        } // end for loop::idxElemCaret

    }// end of forward mapping;


    //doing reverse mapping;
    for(map< Vector3D< size_t >, vector< int > >::iterator it=idxFS.begin(); it!=idxFS.end(); it++)
    {
        const Vector3D<size_t>& key = it->first;
        vector<int>& elemsFreeSurfer = it->second;
        vector<int>& elemsCaret = idxCaret[key];

        if(0==elemsCaret.size()) // counterpart also has elems;
        {
            //search neighbors??? do whole brain , since there are very few;
            WARNING_message("cornered (%d) vertices, will search for whole brain ",elemsFreeSurfer.size());
            elemsCaret.resize(surf1.GetNumOfPoints());
            for ( int idx = 0; idx < surf1.GetNumOfPoints() ; ++idx )
            {
                elemsCaret[idx] = idx;
            } // end for loop::idx

        }

        for ( int idxElemFreeSurfer = 0; idxElemFreeSurfer < elemsFreeSurfer.size() ; ++idxElemFreeSurfer )
        {
            arma::fvec allDistance = arma::randu<arma::fvec>(elemsCaret.size());
            const VectorType& coordFS = surf2.GetPointCoords(elemsFreeSurfer[idxElemFreeSurfer]);
            for ( int idxElemCaret = 0; idxElemCaret < elemsCaret.size() ; ++idxElemCaret )
            {
                const VectorType& coordCaret = surf1.GetPointCoords(elemsCaret[idxElemCaret]);
                VectorType dis = coordFS-coordCaret;
                allDistance[idxElemCaret] = dis.Norm();
            } // end for loop::idxElemFS
            arma::uword minPos=0;
            float minDis = allDistance.min(minPos);
// 	    mapIDCaret2FS[elemsFreeSurfer[idxElemFreeSurfer]] = elemsCaret[minPos];
            mapID2to1[elemsFreeSurfer[idxElemFreeSurfer]] =elemsCaret[minPos];
        } // end for loop::idxElemCaret

    }// end of forward mapping;



}
void CFineParcellatorHCP::SaveIDMapsHFsLv6Native(string anatomy)
{
  if("lh"!=anatomy && "rh"!=anatomy)
    ERROR_exit("unrecognized anatomy: %s, only lh or rh is valid", anatomy.c_str()); 
  
  arma::uvec sub2tmpl(reinterpret_cast<unsigned int*> (&(ndata->mapIDFsLR2FsLv6[0])),ndata->mapIDFsLR2FsLv6.size(),false);
  arma::uvec tmpl2sub(reinterpret_cast<unsigned int*> (&(ndata->mapIDFsLv6toFsLR[0])),ndata->mapIDFsLv6toFsLR.size(),false);
  sub2tmpl.save("subj."+anatomy+".matchAt.tmpl");
  tmpl2sub.save("tmpl."+anatomy+".matchAt.subj");   
}




























};
