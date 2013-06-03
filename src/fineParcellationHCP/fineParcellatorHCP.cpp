/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012.
 *	License : GPL v3
 *
 */

#include "fineParcellatorHCP.h"
#include "string"
#include "triSurface.h"
#include <vector>
#include "ColorSchemeRGB.h"
#include "newimageall.h"
#include <set>
#include "kmlGPMetis.h"
#include "newmat.h"
#include "kmPCA.h"
#include "cfloat"
#include "cmath"
#include "tictoc.hpp"
#include "kaimingCommon.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_rng.h"
#include "omp.h"
#include "armadillo"
#include "gsl/gsl_randist.h"

using namespace KML;
using namespace std;
using namespace NEWIMAGE;
using namespace arma;



namespace KML
{



#define KMLDebug if(DEBUG) std::cout<<__FUNCTION__


class CPrivacyFineParcellatorHCP
{

    friend class CFineParcellatorHCP;

private:

    CPrivacyFineParcellatorHCP ( string surfName ) : surf ( surfName.c_str() ), meanSimilarity ( 0.3 ) ,maxPos4Inter ( 1 ), usePCA4BOLDs(false), randomParcel(false),TR(1.f),numSolutions(50),numTimePoints(0),numPartitions(0){ numFaces=surf.GetNumOfCells(); numVertices=surf.GetNumOfPoints();};
    ~CPrivacyFineParcellatorHCP() {};

    CTriSurface surf;


    //all partions;
    multiset<CPartition> tobePartitions;
    vector<multiset<CPartition> > tmpTobePartitions;
    int maxPos4Inter;
    vector<CPartition> donePartitions;

    //parcels neighbourhood;
    vector<set<int> > nbrParcels;
    vector<vector<int> > nbrCorrelsPcl;

    //partition neighbourhood;
    vector<set<int> > nbrPartitions;
    vector<vector<int> > nbrCorrelsPtn;


    //keep records of all partition lables;
    vector<int> partionLabel4Vertex;


    //the areas of partitions;
    vector<float> allAreas;

    //the cortex label;
    vector<bool> isCortex;
    
    vector<int> allROIVertices;

    //whether pca or average for partition signal;
    bool usePCA4BOLDs;

    //whether or not random parcellation; 
    bool randomParcel; 

    //keep all correlations; 10*abs(correl(i,j));
    vector<vector<int> > allParts;

    vector<float> allRatios;

    CColorSchemeRGB myColorScheme;

    int numVertices;
    int numFaces;
    int numTimePoints;
    int numPartitions;
    int maxNbrs4Parcels;
    int maxNbrs4Partitions;

    float TR;
    float balance;// control via metis;
    float minCorrelThresh;
    int numSolutions;
    float meanSimilarity;

    string mapCaret2FsRegSphere;
    string baseName;
    string fMRIName;

};


CFineParcellatorHCP::CFineParcellatorHCP ( const KML::CFineParcellatorHCP& )
{

}

void CFineParcellatorHCP::PlotPrincipleVector ( std::vector< float >& src )
{
    KML::PlotLine<> ( src );
}

 
void CFineParcellatorHCP::BuildPartitionsAfterNCut()
{

    //build partitions;
    data->donePartitions.clear();
    data->donePartitions.resize ( data->allParts.size() );
    for ( int idxCluster = 0; idxCluster < data->allParts.size() ; ++idxCluster ) {
        data->donePartitions[idxCluster].parcelIds= data->allParts[idxCluster];
        data->donePartitions[idxCluster].numParcels=data->donePartitions[idxCluster].parcelIds.size();
        data->donePartitions[idxCluster].partitionLabel=idxCluster;
        data->donePartitions[idxCluster].numTimePoints=data->numTimePoints;
        this->DoPCA4Partition ( data->donePartitions[idxCluster] );
    } // end for loop::idxCluster

    this->CalculateAreas4Partitions();
//
    this->BuildPartitionNbr();
//
    this->FindPartitionCenters();
//
    //assign color to each partition;
    for ( int idxCluster = 0; idxCluster < data->allParts.size() ; ++idxCluster ) {
        data->donePartitions[idxCluster].color= data->surf.GetPointColor ( data->donePartitions[idxCluster].centerVertexID );
    } // end for loop::idxCluster
};

 
 

void CFineParcellatorHCP::SetROICortexVertices ( string fileName )
{
    data->isCortex.clear();
    data->isCortex.resize ( data->numVertices,false );
    fstream instrm;
    KML::OpenReadStrmAscii ( instrm,fileName );
    string tmpString;
    int numLines=0;
    getline ( instrm,tmpString ); //comment;
    instrm>>numLines;
    cout<<"number of cortex vertices: "<<numLines<<endl;
    for ( int idx = 0; idx < numLines ; ++idx ) {
        int id;
        float x,y,z,alpha;
        instrm>>id>>x>>y>>z>>alpha;
        data->isCortex[id] = true;
    } // end for loop::idx

    instrm.close();
}
void CFineParcellatorHCP::SetRandomParcellation(bool value)
{
  data->randomParcel = value; 
  if(value == true)
  {
   cout<<"Warning: will do random parcellation for the cortex. "<<endl; 
  }

}


CFineParcellatorHCP::CFineParcellatorHCP ( string surfName, string fmriName, string fsLabels )
{
    DEBUG_VAR ( data->baseName );
    data = new CPrivacyFineParcellatorHCP ( surfName, fsLabels );
    data->surf.BuildNeighbors();
    volume4D<float> srcVol;
    NEWIMAGE::read_volume4D ( srcVol, fmriName );

    data->fMRIName = fmriName;
    
    data->randomParcel=false; 
    data->TR = srcVol.TR();
    data->badMapThresh = 10;
    data->numVertices = data->surf.GetNumOfPoints();
    data->numFaces = data->surf.GetNumOfCells();
    data->numTimePoints = srcVol.tsize();
    data->bold4vertices.resize ( data->numVertices, vector<float> ( data->numTimePoints,0 ) );
    data->parcelLabel4vertex.resize ( data->numVertices,-1 );
    data->balance=1000;
    data->maxNbrs4Parcels = 0;
    data->maxNbrs4Partitions=0;
    data->simType = 0;
    data->windLength = 40;
    data->isCortex.resize ( data->numVertices,true );
    data->minCorrelThresh=0.35;
    data->usePCA4BOLDs=false;
    data->numSolutions=100;


}

void CFineParcellatorHCP::SetNumsMetisSolutions(int nSolutions)
{
    data->numSolutions = nSolutions;

}


void CFineParcellatorHCP::SetUsePCA4Bold ( bool usePCA )
{
    data->usePCA4BOLDs= usePCA;
    if(data->usePCA4BOLDs)
        cout<<"use pca for partition bolds"<<endl;
    else
        cout<<"use mean for partion bolds"<<endl;
}

void CFineParcellatorHCP::SetMinCorrelThresh ( float minThresh )
{
    KML::CheckValueRange<> ( minThresh,0.f,1.f );
    data->minCorrelThresh = minThresh;
}

float CFineParcellatorHCP::GetMinCorrelThresh ( void ) const
{
    return data->minCorrelThresh;
}


void CFineParcellatorHCP::SetCaretSurf2FSReg ( string name )
{
    data->mapCaret2FsRegSphere = name;
    //TODO::add read file; 

}
int CFineParcellatorHCP::SetBaseName ( string name )
{
    data->baseName = name;

}
  
float CFineParcellatorHCP::GetBalance ( void )
{
    return data->balance;
}
void CFineParcellatorHCP::SetBalance ( float balance )
{
    data->balance = balance;
}
  

void CFineParcellatorHCP::FindPartitionCenters()
{
    // only do this after partiion is done;
    for ( int idxPart = 0; idxPart < data->donePartitions.size() ; ++idxPart ) {
        VectorType meanCoord;
        int count = data->donePartitions.operator[] ( idxPart ).parcelIds.size();

        for ( int idxPartVertices = 0; idxPartVertices < data->donePartitions.operator[] ( idxPart ).parcelIds.size() ; ++idxPartVertices ) {
            int vid = data->donePartitions.operator[] ( idxPart ).parcelIds[idxPartVertices];
            meanCoord += data->surf.GetPointCoords ( vid );
        } // end for loop::idxPcls

        meanCoord/= count;

        //find the closest vertex;
        int closeVertexId = -1;
        float minDistance = 10000;

        for ( int idxPartVertices = 0; idxPartVertices < data->donePartitions.operator[] ( idxPart ).parcelIds.size() ; ++idxPartVertices ) {
            int vid = data->donePartitions.operator[] ( idxPart ).parcelIds[idxPartVertices];
            const VectorType& theCoord = data->surf.GetPointCoords ( vid );
            VectorType disCoord = theCoord - meanCoord;
            if ( disCoord.NormSquare() < minDistance ) {
                closeVertexId  = vid;
                minDistance = disCoord.NormSquare();
            }
        } // end for loop::idxPcls

        if ( -1==closeVertexId ) {
            cout<<"error while finding centerVertexID"<<endl;
            exit ( 1 );
        } else {
            data->donePartitions[idxPart].centerVertexID = closeVertexId;
        }
    } // end for loop::idxPart


//cout for checking;
//     for ( int idx = 0; idx < data->donePartitions.size() ; ++idx )
//     {
//         cout<<idx<<" "<<data->donePartitions[idx].centerVertexID<<endl;
//     } // end for loop::idx


}

int CFineParcellatorHCP::GetMaxNbrPartitions ( void )
{

    if ( data->maxNbrs4Partitions ) {
        return data->maxNbrs4Partitions;
    }

    int maxNum = 0;
    int thePartId = -1;
    for ( int idxPartitions = 0 ; idxPartitions < data->nbrPartitions.size(); ++idxPartitions ) {
        if ( data->nbrPartitions[idxPartitions].size() > maxNum ) {
            maxNum= data->nbrPartitions[idxPartitions].size();
            thePartId = idxPartitions;
        }
    } // end of for::idxPartitions

    cout<<"maximum number of neighbour partitions: "<<maxNum<<endl;
//     this->PrintPartition(data->donePartitions[thePartId]);
    data->maxNbrs4Partitions = maxNum+1;
    return data->maxNbrs4Partitions;

}

 
fcol& CFineParcellatorHCP::GetVertexTimeSeries ( int vertex )
{
    return data->parcels.operator[] ( parcelIdx ).timeSeries;
}
float CFineParcellatorHCP::GetTR ( void )
{
    return data->TR;
}
int CFineParcellatorHCP::GetNumFaces ( void )
{
    return data->numFaces;
} 
int CFineParcellatorHCP::GetNumVertices ( void )
{
    return data->numVertices;
}

void CFineParcellatorHCP::SetTR ( float tr )
{
    data->TR = tr;
}
  
CFineParcellatorHCP::~CFineParcellatorHCP()
{
    delete data;
}

void CFineParcellatorHCP::SavePartitionsOnSurface ( string fileName )
{

    data->surf.NeedPointsColor ( );
    data->surf.NeedPointsLabel ( );
    for ( int idx = 0 ; idx <data->numVertices; ++idx ) {
        int parcelLabel = data->partionLabel4Vertex[ idx];
        data->surf.SetPointColor ( idx,data->donePartitions[parcelLabel].color );
        data->surf.SetPointLabel ( idx,data->partionLabel4Vertex[idx] );
    } // end of for::idx

    data->surf.SaveAs ( fileName );

}

void CFineParcellatorHCP::PrintPartition ( const CPartition& currentPartition )
{
//     const_cast<CPartition&>(currentPartition)>>cout;

    cout<<"partition id: "<<currentPartition.partitionLabel<<endl<<"number of parcels: "<<currentPartition.numParcels<<endl;
//         <<"neighbour parts: "<<currentPartition.nbrPartitions<<endl<<"principalVector: "<<currentPartition.principalVector<<endl
//         <<"center vertex id: "<<currentPartition.centerVertexID<<endl;

}
   
void CFineParcellatorHCP::CheckZeroPrincipalVectors4Parts ( void )
{
    int numZeroVectors = 0;

    for ( int idx = 0; idx < data->donePartitions.size() ; ++idx ) {
        bool isAllZeros = true;

        for ( int t = 1; t < data->numTimePoints ; ++t ) {
            if ( data->donePartitions[idx].principalVector[t] != 0 ) {
                isAllZeros = false;
                break;
            }
        } // end for loop::t

        if ( isAllZeros ) {
            ++numZeroVectors;
            cout<<"##############All Zero Principle Vectors:"<<idx<<endl;
            this->PrintPartition ( data->donePartitions[idx] );
        }
    } // end for loop::idx

    if ( numZeroVectors==0 ) {
        cout<<"################## No Zero Vectors Found ###############" <<endl;
    } else {
        cout<<"##################"<< numZeroVectors<<" Zero Vectors Found!############"<<endl;
    }
}

void CFineParcellatorHCP::PartitioningWholeBrainWithMetis ( int nParts )
{

    data->baseName+=".";
    data->baseName+=nParts;

    //indexing surface;
    this->BuildInitParcels ();

    //do nc;
    this->DoPartitioning ( nParts );

    //build partitions;
    this->BuildPartitionsAfterNCut();
// //
    if ( data->mapCaret2FsRegSphere.size() ) {
        this->SavePartitionMetaInforSphere ( data->mapCaret2FsRegSphere,data->baseName );
    }
// //
    this->PrintAllPartitionInfo ( false,data->baseName+".partInfo" );

    cout<<data->allRatios<<endl;
    cout<<KML::StatsMean(data->allRatios)<<endl;
    cout<<KML::StatsStdev(data->allRatios)<<endl;

}

float CFineParcellatorHCP::ComputeNCut(graph_t *graph, int *where, int npart)
{
    int i, j, cm, nvtxs;
    idx_t *ncut, *degree, *xadj, *adjncy;
    float result;
    idx_t * adjwgt;

    ncut = ismalloc(npart, 0, "ComputeNCut: ncut");
    degree = ismalloc(npart, 0, "ComputeNCut: degree");
    nvtxs = graph->nvtxs;
    xadj = graph->xadj;
    adjncy = graph->adjncy;
    adjwgt = graph->adjwgt;

    if (graph->adjwgt == NULL) {
        for (i=0; i<nvtxs; i++) {
            cm = where[i];
            for (j=xadj[i]; j<xadj[i+1]; j++) {
                degree[cm] ++;
                if (cm != where[adjncy[j]])
                    ncut[cm] ++;
            }
        }
    }
    else {
        for (i=0; i<nvtxs; i++) {
            cm = where[i];
            for (j=xadj[i]; j<xadj[i+1]; j++) {
                degree[cm] += adjwgt[j];
                if (cm != where[adjncy[j]])
                    ncut[cm] += adjwgt[j];
            }
        }
    }
    int empty = 0;
    result =0;
    for (i=0; i<npart; i++) {
        if (degree[i] == 0)
            empty++;
        if (degree[i] >0)
            result +=  ncut[i] *1.0/ degree[i];
    }
    //printf("Empty clusters: %d\n", empty);
    free(ncut);
    free(degree);
    return result+empty;
}


/*************************************************************************
* This function computes the cut given the graph and a where vector
**************************************************************************/
int CFineParcellatorHCP::ComputeGraphCut(graph_t* graph, int* where)
{
    int i, j, cut;

    if (graph->adjwgt == NULL) {
        for (cut=0, i=0; i<graph->nvtxs; i++) {
            for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
                if (where[i] != where[graph->adjncy[j]])
                    cut++;
        }
    }
    else {
        for (cut=0, i=0; i<graph->nvtxs; i++) {
            for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
                if (where[i] != where[graph->adjncy[j]])
                    cut += graph->adjwgt[j];
        }
    }

    return cut/2;
}


void CFineParcellatorHCP::DoPartitioning ( int numCuts )
{

    graph_t* currentGraph = CreateGraph();
    this->PrepareGraph4Metis ( currentGraph );

    params_t* params = new params_t;
    KmlInitParams(params,numCuts);
    params->niter=10;
    params->ncuts=data->numSolutions;
    params->ufactor=data->balance; // that allows the maximal size is twice of average size;

    vector<int> newPartitionsInfo;

    KmlGraphPartitionMetis ( currentGraph,params, newPartitionsInfo, true);

    vector<vector<int> >& allParts=data->allParts;
    allParts.resize ( numCuts );

    data->donePartitions.clear();
    for ( int idx = 0; idx < newPartitionsInfo.size() ; ++idx ) {
        // newPartitionsInfo contains node label;
        int label = newPartitionsInfo[idx];
        int pid = data->allROIVertices[idx];
        allParts[label].push_back ( pid );
    } // end for loop::idx

    data->surf.NeedPointsLabel();
    data->surf.NeedPointsColor();
    for ( int idx = 0; idx < newPartitionsInfo.size() ; ++idx ) {
        int label = newPartitionsInfo[idx];
        int pid = data->allROIVertices[idx];
        data->surf.SetPointLabel ( pid,label );
    } // end for loop::idx

    //check the partitions that has problem;
//     this->CheckScattersAndFix();
    data->surf.ColorizeSurfWLabels();


    cout<<"normalized cut: "<< this->ComputeNCut(currentGraph, &newPartitionsInfo[0],numCuts)<<endl;
    cout<<"graph cut: "<< this->ComputeGraphCut(currentGraph, &newPartitionsInfo[0])<<endl;

    delete params ;
    FreeGraph(&currentGraph);
}

void CFineParcellatorHCP::SavePartitions ( bool toSaveVtk, int pid )
{
    //check pcc connectivity pattern if pid set;
    if ( -1!=pid ) {
        int pccLabel = data->surf.GetPointLabel ( pid );
        vector<float> allCorrels ( data->donePartitions.size(),0 );
        for ( int idxLabel = 0; idxLabel < data->donePartitions.size() ; ++idxLabel ) {
            allCorrels[idxLabel] = this->GetSimilarityTwoTimeSeries ( data->donePartitions[idxLabel].principalVector,data->donePartitions[pccLabel].principalVector );
            allCorrels[idxLabel]= ( allCorrels[idxLabel] );
        } // end for loop::idxLabel

        vector<float> attriCorrel ( data->numVertices,0 );
        for ( int idxV = 0; idxV < data->numVertices ; ++idxV ) {
            int label= data->surf.GetPointLabel ( idxV );
            if ( -1 != label ) {
                attriCorrel[idxV] = allCorrels[label];
            }
        } // end for loop::idxV
        string attiName = "correl.";
        attiName+=pid;
        data->surf.AddPointDataScalar ( attriCorrel,attiName );
    }
    
    if(toSaveVtk)
    {
      string fileName = data->baseName+".vtk";
      cout<<"write partitions to vtk file: "<< fileName<<endl;
      data->surf.SaveAs ( fileName );
    }

}
 
void CFineParcellatorHCP::RebuildPartitionsAndNbr()
{

    //build new nbrs and partitions;
    cout<<"num of original partitions:"<< data->donePartitions.size() <<endl;
    vector<CPartition> newPartitions;
    for ( int idx = 0; idx < data->nbrPartitions.size() ; ++idx ) {
        if ( data->donePartitions[idx].numParcels ) {
            newPartitions.push_back ( data->donePartitions[idx] );
        }
        if ( data->donePartitions[idx].numParcels && data->nbrPartitions[idx].size() ==0 ) {
            cout<<"Error here: "<<idx<<endl;
        }

    } // end for loop::idx

    data->donePartitions = newPartitions;
    data->numPartitions = data->donePartitions.size();
    cout<<"num of new partitions:"<< data->donePartitions.size() <<endl;

    this->BuildPartitionNbr();

}

void CFineParcellatorHCP::DoPCA4Partition ( CPartition& currentPartition )
{

    set<int> parcelIds;
    for ( int idx = 0; idx < currentPartition.numParcels; ++idx ) {
        int pid = currentPartition.parcelIds[idx];
        int parcelId = data->parcelLabel4vertex[pid];
        parcelIds.insert ( parcelId );
    } // end for loop::idx_t
    std::set< int >::iterator it=parcelIds.begin();

    //will do pca anyway;

    currentPartition.principalVector.clear();
    currentPartition.principalVector.resize ( data->numTimePoints,0);

    CPCA pca;
    vector<vector<float> > rawdata, basis;
    vector<float> eigValue;
    while ( it++!=parcelIds.end() ) {
        rawdata.push_back ( data->parcels[*it].timeSeries );
    }
    pca.SetData ( rawdata );
    pca.DoPCA ( basis,eigValue );

    float sumEvalues = KML::StatsSum ( eigValue );
    float ratio=0;
    if ( 0!=sumEvalues ) {
        currentPartition.ratio = eigValue[0]/ ( sumEvalues+FLT_MIN );
    }
    else {        //for sigular pca;
        cout<<"Warnig: singular pca for time series"<<endl;
        currentPartition.ratio=1;
    }

    data->allRatios.push_back(currentPartition.ratio);

    if ( data->usePCA4BOLDs ) {

        for ( int idx = 0; idx < data->numTimePoints ; ++idx ) {
            currentPartition.principalVector[idx]= basis[idx][0];
        } // end for loop::idx

    } else {

        arma::fvec averageTimeSeries ( data->numTimePoints );
        averageTimeSeries.fill ( 0 );
        it=parcelIds.begin();
        while(it!=parcelIds.end())
        {
            vector<float>& timeSeries= data->parcels[*it].timeSeries;
            arma::fvec armaTS ( & ( timeSeries[0] ),timeSeries.size() );
            averageTimeSeries+=armaTS;
            it++;
        }

        averageTimeSeries/= ( float ( parcelIds.size() ) );
        currentPartition.principalVector.resize ( data->numTimePoints,0 );
        for ( int idx = 0; idx < data->numTimePoints ; ++idx ) {
            currentPartition.principalVector[idx]=averageTimeSeries[idx];
        } // end for loop::idx
    }


    //output all time series;
//     {
//       vector<vector<float> > rawdata;
//       it=parcelIds.begin();
//       while ( it++!=parcelIds.end() ) {
//             rawdata.push_back ( data->parcels[*it].timeSeries );
//       }
//       rawdata.push_back(currentPartition.principalVector);
// 
//       string fileName("partition.");
//       fileName+=currentPartition.partitionLabel;
//       KML::SaveFloatMatrix(fileName,rawdata);
//     }
}
 
void CFineParcellatorHCP::CheckPartitionTopology()
{
    //for each label;
    //find all vertices, and marks in a vector;
    //start with a vertex, expand to max range with the same label, find all the vertices,
    //compare the two groups of vertices if not equal, then has topology issue;
    //for the minority, use its neighborhood lable information to relabel the vertex;
//   data->surf.BuildNeighbors();

    vector<int>& allLabels = data->surf.GetAllPointLabel();
//     cout<<allLabels<<endl;

    for ( int idxLabel = 0; idxLabel < data->donePartitions.size() ; ++idxLabel ) {
        vector<int> allVWLabel;
        vector<int> maxHomoVWLabel;

        for ( int idxV = 0; idxV < allLabels.size() ; ++idxV ) {
            if ( idxLabel == allLabels[idxV] ) {
                allVWLabel.push_back ( idxV );
            }
        } // end for loop::idxV


        if ( allVWLabel.size() ==0 ) {
            cout<<"Warning: Non existence partition:  "<< idxLabel<<endl;
            continue;
        }

        list<int> groupInfos;
        map<int, list<int> > nbrInfos;

        data->surf.Propagate2All ( allVWLabel[0],groupInfos,nbrInfos,allLabels );

        if ( groupInfos.size() != allVWLabel.size() ) {
            cout<<"Partition has topology issue: "<< idxLabel<<" "<<allVWLabel.size() <<" "<<groupInfos.size() <<endl;;
            //suppose there are only two groups, determine which group is smaller, then will kill the smaller one;
            if ( groupInfos.size() < 0.5*allVWLabel.size() ) {

            } else {
                //fine one that not in the first group, and kill use its neighbour information to kill it;

            }

        }




    } // end for loop::idxLabel



}


void CFineParcellatorHCP::BuildPartitionNbr ( void )
{

    cout<<"BuildPartitionNbr"<<endl;
    data->numPartitions = data->donePartitions.size();
    data->partionLabel4Parcels.clear();
    data->partionLabel4Vertex.clear();
    data->partionLabel4Parcels.resize ( data->numParcels,-1 );
    data->partionLabel4Vertex.resize ( data->numVertices,-1 );

    #pragma omp parrel for
    // generate partition label;
    for ( int idxPartitions = 0; idxPartitions < data->donePartitions.size() ; ++idxPartitions ) {

        //for each parcel, assign partition label ;
        for ( int idxPartVertex = 0; idxPartVertex < data->donePartitions.operator[] ( idxPartitions ).parcelIds.size() ; ++idxPartVertex ) {
            int pid = data->donePartitions[idxPartitions].parcelIds[idxPartVertex];
            data->partionLabel4Vertex[pid] = idxPartitions;
        } // end for loop::idxParcels
    } // end for loop::idxPartitions

    //build partion neighborhood;
    data->nbrPartitions.clear();
    data->nbrPartitions.resize ( data->numPartitions );

    for ( int idxV = 0; idxV < data->numVertices ; ++idxV ) {
        int label1 = data->partionLabel4Vertex[idxV];
        if ( -1!=label1 ) { // -1 for non interested rois, i.e., medial wall;
            list<int> nbrs;
            data->surf.GetNeighboringPointsOfPoint ( idxV,nbrs );
            std::list< int >::iterator it=nbrs.begin();
            while ( it!=nbrs.end() ) {
                int label2 = data->partionLabel4Vertex[*it];
                if ( -1!=label2 && label2!=label1 ) {
                    data->nbrPartitions[label2].insert ( label1 );
                    data->nbrPartitions[label1].insert ( label2 );
                }
                it++;
            }
        }

    } // end for loop::idxV


//     cout<<"output nbr information"<<endl;
//     for ( int idxL = 0; idxL < data->donePartitions.size() ; ++idxL )
//     {
//         cout<<"label/vertex"<<idxL<<"/"<<data->donePartitions[idxL].parcelIds[0]<<"  :::  ";
//         std::set< int >&  nbrs = data->nbrPartitions[idxL];
//         std::set< int >::iterator it=nbrs.begin();
//         while(it!=nbrs.end())
//         {
//             int label = *it;
//             cout<<data->donePartitions[label].parcelIds[0]<<" ";
//             it++;
//         }
//         cout<<endl;
//     } // end for loop::idxL

}

void CFineParcellatorHCP::GenerateColors4Partitions()
{


    //generate colors for partitions;
    int maxColorNums = this->GetMaxNbrPartitions();
    data->myColorScheme.SetColorNums ( maxColorNums );
    data->donePartitions.operator[] ( 0 ).color = data->myColorScheme.GetColorByIndex ( 0 );

    for ( int idxPartition = 1; idxPartition < data->donePartitions.size() ; ++idxPartition ) {
        list<RGBTYPE> usedColors;
        for ( set<int>::iterator itNbrPartions = data->nbrPartitions.operator[] ( idxPartition ).begin(); itNbrPartions != data->nbrPartitions[idxPartition].end(); itNbrPartions++ ) {
            if ( *itNbrPartions < idxPartition ) {
                usedColors.push_back ( data->donePartitions[*itNbrPartions].color );
            }
        }

        data->donePartitions[idxPartition].color = data->myColorScheme.GetRandomColorWithExclusionS ( usedColors );
    } // end for loop::idxPartition

}

void CFineParcellatorHCP::PrepareGraph4Metis ( graph_t* graph, const KML::CPartition& currentPartition )
{

    graph->nedges=0;

    vector<bool> selectedParcels ( data->numParcels, false );
    vector<int> parcelNodeID ( data->numParcels,0 );

    for ( int idx = 0; idx < currentPartition.numParcels ; ++idx ) {
        int parcelIdx =currentPartition.parcelIds[idx];
        selectedParcels[ parcelIdx ] = true;
        parcelNodeID[parcelIdx] = idx;
    } // end for loop::idx


    for ( int idx = 0; idx < currentPartition.numParcels ; ++idx ) {
        int currendParcel = currentPartition.parcelIds[idx];
        for ( set<int>::iterator itNbrs = data->nbrParcels [currendParcel].begin(); itNbrs != data->nbrParcels[currendParcel].end(); itNbrs++ ) {
            if ( selectedParcels[*itNbrs] ) {
// 	cout<<parcelNodeID[*itNbrs]<<" ";
                ++ ( graph->nedges );
            }
        }
//     cout<<endl;
    } // end for loop::idx

//   cout<<parcelNodeID<<endl;
    graph->nvtxs= currentPartition.parcelIds.size();
//     graph->ncon=1;
    graph->xadj = ismalloc ( graph->nvtxs+1, 0, const_cast<char*> ( string ( "ReadGraph: xadj" ).c_str() ) );
    graph->adjncy = imalloc ( graph->nedges, const_cast<char*> ( string ( "ReadGraph: adjncy" ).c_str() ) );
    graph->adjwgt =  ismalloc ( graph->nedges, 1,const_cast<char*> ( string ( "ReadGraph: adjwgt" ).c_str() ) );




//   xadj   = graph->xadj   = ismalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
//   adjncy = graph->adjncy = imalloc(graph->nedges, "ReadGraph: adjncy");
//   vwgt   = graph->vwgt   = ismalloc(ncon*graph->nvtxs, 1, "ReadGraph: vwgt");
//   adjwgt = graph->adjwgt = ismalloc(graph->nedges, 1, "ReadGraph: adjwgt");
//   vsize  = graph->vsize  = ismalloc(graph->nvtxs, 1, "ReadGraph: vsize");


    //set values;

    int idxEdge=0;

    graph->xadj[0]=0;
    for ( int idxNode = 0; idxNode < currentPartition.numParcels ; ++idxNode ) {
        int currendParcel = currentPartition.parcelIds[idxNode];
        int currentParcelNodeIdx = idxNode;

        int idxNbrParcelInALLNeighbors = 0;  // keep record the pos in nbr parcesl;
        for ( set<int>::iterator itNbrs = data->nbrParcels [currendParcel].begin(); itNbrs != data->nbrParcels[currendParcel].end(); itNbrs++ ) {
            int nbrParcel =*itNbrs;
            if ( selectedParcels[nbrParcel] ) {
// 	cout<<data->nbrCorrels[currendParcel][idxNbrParcelInALLNeighbors]<<" ";
                int nbrParcelNodeIdx = parcelNodeID[nbrParcel];
                graph->adjncy[idxEdge]= nbrParcelNodeIdx;
                graph->adjwgt[idxEdge]= data->nbrCorrelsPcl[currendParcel][idxNbrParcelInALLNeighbors];
                ++idxEdge;
            }

            ++idxNbrParcelInALLNeighbors;
        }
//     cout<<endl;
        graph->xadj[idxNode+1] = idxEdge;
    } // end for loop::idx

}

void CFineParcellatorHCP::PrepareGraph4Metis ( graph_t* graph )
{
  
  /**** Start of section::for random edge;  ****/
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
          

    graph->nedges=0;
    graph->nvtxs=0;

    vector<int>& allROIVertices= data->allROIVertices;
    allROIVertices.clear();
    vector<int> pid2nodeIdx ( data->numVertices,-1 ); //map pid to nid;

    for ( int idx = 0; idx < data->numVertices ; ++idx ) {
        if ( data->isCortex[idx] ) {
            allROIVertices.push_back ( idx );
            list<int> nbrs;
            data->surf.GetNeighboringPointsOfPoint ( idx,nbrs );
            for ( list<int>::iterator itNbrs = nbrs.begin(); itNbrs != nbrs.end(); itNbrs++ ) {
                if ( data->isCortex[*itNbrs] ) {
                    ++ ( graph->nedges );
                }
            }
        }
    } // end for loop::idx

    for ( int idxROIV = 0; idxROIV < allROIVertices.size() ; ++idxROIV ) {
        pid2nodeIdx[allROIVertices[idxROIV]]=idxROIV;
    } // end for loop::idxROIV

//   cout<<parcelNodeID<<endl;
//     graph->ncon=1;
    graph->nvtxs=allROIVertices.size();
    graph->ncon=1;
    graph->nedges*=2;
    graph->xadj = ismalloc ( graph->nvtxs+1, 0, const_cast<char*> ( string ( "ReadGraph: xadj" ).c_str() ) );
    graph->adjncy = imalloc ( graph->nedges, const_cast<char*> ( string ( "ReadGraph: adjncy" ).c_str() ) );
    graph->adjwgt =  ismalloc ( graph->nedges, 1,const_cast<char*> ( string ( "ReadGraph: adjwgt" ).c_str() ) );
    graph->vwgt   = ismalloc(graph->nvtxs, 1, "ReadGraph: vwgt");


    //set values;

    int idxEdge=0;
    graph->xadj[0]=0;
    for ( int idxNode = 0; idxNode < allROIVertices.size() ; ++idxNode ) {



        int pointId = allROIVertices[idxNode];
        list<int> nbrs;
        data->surf.GetNeighboringPointsOfPoint ( pointId,nbrs );
        for ( list<int>::iterator itNbrs = nbrs.begin(); itNbrs != nbrs.end(); itNbrs++ ) {
            //determine weights;
            int nbrPid = *itNbrs;
            int nbrNid = pid2nodeIdx[nbrPid];
            int weight = 1;

            //if both are not wall; then regular;
            //else if both on wall; then high;
            //esle boarder; then 0;
            if ( data->isCortex[pointId] && data->isCortex[nbrPid] ) {

	        graph->adjncy[idxEdge]=nbrNid;
		if(data->randomParcel) 
		{
		  double k = abs(gsl_ran_gaussian (r, 0.2))+0.2;
		  graph->adjwgt[idxEdge]= k*100;
		}
		else
		{         
		  int parcel1 = data->parcelLabel4vertex[pointId];
		  int parcel2 = data->parcelLabel4vertex[nbrPid];
		  float correl = this->GetSimilarityTwoTimeSeries ( data->parcels[parcel1].timeSeries, data->parcels[parcel2].timeSeries );
		  if(correl>0)
		    weight = int ( abs ( correl ) *100 );
		  graph->adjwgt[idxEdge]=weight;
		}		 
                ++idxEdge;
            }

        }
        graph->xadj[idxNode+1] = idxEdge;

    }
    
    
    gsl_rng_free (r);

}


void CFineParcellatorHCP::SaveGraph ( graph_t* graph, string fileName )
{

    WriteGraph(graph,const_cast<char*> (fileName.c_str()));
 
}


void CFineParcellatorHCP::PrintAllPartitionInfo ( bool isInteractive, string saveFile )
{

    vector<double> allRatios;
    vector<double> allSizes;
    uvec histRatio,histSizes;
    {
        //set values for histgram;
        for ( int idx = 0; idx < data->numPartitions ; ++idx ) {
            allRatios.push_back ( data->donePartitions[idx].ratio );
            allSizes.push_back ( data->donePartitions[idx].parcelIds.size() );
        } // end for loop::idx
        vec allRatiosArma ( &allRatios[0],allRatios.size(),false );
        vec allSizesArma ( &allSizes[0],allSizes.size(),false );
        histRatio = hist ( allRatiosArma,linspace<vec> ( 0,1,101 ) );
        histSizes = hist ( allSizesArma,linspace<vec> ( 0,2000,101 ) );
    }



    if ( saveFile.size() ) {
        fstream outStrm;
        KML::OpenWriteStrmAscii ( outStrm,saveFile );
        outStrm<< "number of partitions:\t"<< data->donePartitions.size() <<endl;
        outStrm<< "histogram of PCA ratio:"<<endl;
        outStrm<< histRatio.t() <<endl;
        outStrm<< "histogram of partition size(number of vertices):"<<endl;
        outStrm<< histSizes.t() <<endl;
        outStrm<< "areas of partitions: "<<endl;
        outStrm<< data->allAreas<<endl;
        outStrm.close();
    }
}

void CFineParcellatorHCP::SavePartitionMetaInforSphere ( string sphereReg, string fileName )
{
    CTriSurface mySurf ( sphereReg,true );
    float radius = 0; //default freesurfer sphere radius;
    VectorType center;

    //estimate center;
    for ( int idx = 0; idx < mySurf.GetNumOfPoints() ; ++idx ) {
        const VectorType& coord = mySurf.GetPointCoords ( idx );
        center+=coord;
    } // end for loop::idx
    center/=mySurf.GetNumOfPoints();
    cout<<"Center coord: "<< center<<endl;

    //estimate radius;
    for ( int idx = 0; idx < mySurf.GetNumOfPoints() ; ++idx ) {
        const VectorType& coord = mySurf.GetPointCoords ( idx );
        VectorType dis = coord - center;
        float tmpR = dis.Norm();
        radius = tmpR > radius ? tmpR : radius;
    } // end for loop::idx

    cout<<"Radius: "<<radius<<endl;

    //get the r theta, phi; information;
    vector<VectorType> coordsSurfs;
    vector<VectorType> coordsSpheres;

    for ( int idxPart = 0; idxPart < data->donePartitions.size() ; ++idxPart ) {
        const VectorType& coord = data->surf.GetPointCoords ( data->donePartitions[idxPart].centerVertexID );
        coordsSurfs.push_back ( coord );
        coordsSpheres.push_back ( mySurf.GetPointCoords ( data->donePartitions[idxPart].centerVertexID ) );

        VectorType standardCoord = mySurf.GetPointCoords ( data->donePartitions[idxPart].centerVertexID )-center;

        float thitaValue = standardCoord.z/radius;
        if ( thitaValue > 1 || thitaValue < -1 ) {
            cout<<"Warning : Invalid thitaValue :"<< thitaValue<<endl;
            thitaValue = thitaValue > 1 ? 1 : thitaValue;
            thitaValue = thitaValue <-1 ? -1 : thitaValue;
        }

        data->donePartitions[idxPart].sphereCoords.x = radius;
        data->donePartitions[idxPart].sphereCoords.y = acos ( thitaValue );
        data->donePartitions[idxPart].sphereCoords.z = atan2 ( standardCoord.y,standardCoord.x );
//     cout<<data->donePartitions[idxPart].sphereCoords.y*180/3.1415926<<" "<<data->donePartitions[idxPart].sphereCoords.z*180/3.1415926<<endl;
    } // end for loop::idxPart

    string baseName= fileName;
    //write center as vertices;
    // on surface;
    string centerOnSurfName = baseName + ".CentersonSurf.vtk";
    KML::VisualizePointsUsingVerticeByVTK ( coordsSurfs, centerOnSurfName );
    //write center as vertices;
    // on sphere;
    string centerOnSphereName = baseName + ".CentersonSphere.vtk";
    KML::VisualizePointsUsingVerticeByVTK ( coordsSpheres,centerOnSphereName );

    //write metainfo
    fstream outStrm;
    string metaName = baseName+".CenterMetaInfo";

    KML::OpenWriteStrmBinary ( outStrm,metaName );
    for ( int idxPart = 0; idxPart < data->donePartitions.size() ; ++idxPart ) {
        //write nbr partitions;
        data->donePartitions[idxPart].nbrPartitions.clear();
        data->donePartitions[idxPart].numNbrParts=data->nbrPartitions[idxPart].size();

        set<int>::iterator it = data->nbrPartitions[idxPart].begin();
        while ( it != data->nbrPartitions[idxPart].end() ) {
            data->donePartitions[idxPart].nbrPartitions.push_back ( *it );
            it++;
        }
        data->donePartitions[idxPart]>>outStrm;
    } // end for loop::idxPart

    outStrm.close();

//     vector<VectorType> allPoints;
//     set<int> allNbrParts = data->nbrPartitions[0];
//     cout<<allNbrParts<<endl;
//     cout<<data->donePartitions[0].nbrPartitions<<endl;
//     for(std::set< int >::iterator it = allNbrParts.begin(); it!= allNbrParts.end(); it++)
//     {
//       int partitionLabel = *it;
//
//
//     for ( int idx = 0; idx < data->surf.GetNumOfPoints() ; ++idx )
//     {
//       if(partitionLabel == data->surf.GetPointLabel(idx))
// 	allPoints.push_back(data->surf.GetPointCoords(idx));
//     } // end for loop::idx
//
//     }
//
//     KML::VisualizePointsUsingVerticeByVTK(allPoints,"checkNbrOri.vtk");

//see whether read works;
//   vector<CPartition> tmpParts;
//   KML::OpenReadStrmBinary(outStrm, metaName);
//   for ( int idx = 0; idx < data->numPartitions ; ++idx )
//   {
//     CPartition currpart;
//     currpart<<outStrm;
//     if( currpart!= data->donePartitions[idx])
//       cout<<"ERROR! Unmatched: "<<idx<<endl;
//     else
//       cout<<"matched: "<<idx<<endl;
//   } // end for loop::idx
//   outStrm.close();



}

void CFineParcellatorHCP::OutputPartitionSizeHistogram ( string fileName, int minSize, int maxSize )
{

    vector<int> hist ( maxSize-minSize,0 );
    for ( int idxPart = 0; idxPart < data->donePartitions.size() ; ++idxPart ) {
        ++hist[ data->donePartitions[idxPart].numParcels- minSize];
    } // end for loop::idxPart

    KML::SaveStdVectorAsLine<> ( fileName,hist );

}

float CFineParcellatorHCP::GetSimilarityTwoTimeSeries ( std::vector< float >& timeseries1, std::vector< float >& timeseries2 )
{
    if ( 0== data->simType ) {
        float correl = KML::CoorelOfTwoSeries<> ( timeseries1,timeseries2,0,data->numTimePoints ) ;
        correl = abs ( correl ) > data->minCorrelThresh ? correl : 0;
        return correl;

    } else if ( 1 == data->simType ) {

        // get the minimun value of sliding windows;
        vector<float> allSims;
        for ( int start = 0; start < data->numTimePoints-data->windLength ; ++start ) {
            float correl = KML::CoorelOfTwoSeries<> ( timeseries1,timeseries2, start, start+ data->windLength );
            correl = abs ( correl );
            allSims.push_back ( correl );
        } // end for loop::start

        float minValue = 100;
        int minPos = -1;

        KML::findMinValueAndIndex<> ( allSims, minValue, minPos );
        if ( -1==minPos | minPos >= allSims.size() ) {
            cout<<"Error! find min value and index for sliding window."<<endl;
            exit ( 1 );
        }
        return minValue;
    } else {

        cout<<"Error! Unsupported similarity type. "<<endl;
        exit ( 1 );
    }

    return 0;

}



void CFineParcellatorHCP::CalculateInfoCritera ( float& part1, float& aic, float& bic, float& icl )
{
    //estimate distribution variance;
    //the variance sigmaSquare = 1/N*SUM(SUM(Dis(pi-Ck)); pi is a member of cluster/partition Ck; N is the number of parces;
    float sigmaSquare = 0;
    float meanVar = 0;

    vector<float> vMeanDis ( data->donePartitions.size(),0 );
    vector<float> vVarianceDis ( data->donePartitions.size(),0 );

    int numberParameters = data->donePartitions.size() * ( 1/*data->numTimePoints-2*/ ) + 1;
    for ( int idxCluster = 0; idxCluster < data->donePartitions.size() ; ++idxCluster ) {
        for ( int idxParcels = 0; idxParcels < data->donePartitions[idxCluster].numParcels ; ++idxParcels ) {
            float dis = this->GetDistanceHParcelAndPartition ( data->donePartitions[idxCluster].parcelIds[idxParcels], idxCluster );
            meanVar += dis;
//       vMeanDis[idxCluster]+=this->GetSqureDistanceHParcelAndPartition(data->donePartitions[idxCluster].parcelIds[idxParcels], idxCluster);
        } // end for loop::idxParcels
//     vMeanDis[idxCluster]/=data->donePartitions[idxCluster].numParcels;
    } // end for loop::idxCluster
    meanVar/=data->numParcels;


    vector<float> allRatios ( data->donePartitions.size(),0 );
    for ( int idxCluster = 0; idxCluster < data->donePartitions.size() ; ++idxCluster ) {
        //for each cluster, estimate the variance;
        vector<float> allDisWithin, allDisNbr;
        CPartition& partition = data->donePartitions[idxCluster];
        //measure within partition;
        for ( int idxParcel = 0; idxParcel < data->donePartitions[idxCluster].numParcels ; ++idxParcel ) {
            int parcelID = data->donePartitions[idxCluster].parcelIds[idxParcel];
            float disSample = this->GetDistanceHParcelAndPartition ( data->donePartitions[idxCluster].parcelIds[idxParcel], idxCluster );
            allDisWithin.push_back ( disSample );
        } // end for loop::idxParcel//the variance sigmaSquare = 1/N*SUM(SUM(Dis(pi-Ck)); pi is a member of cluster/partition Ck; N is the number of parces;

        const std::set< int >& nbrPartitions = data->nbrPartitions[idxCluster];

        for ( std::set< int >::iterator it = nbrPartitions.begin(); it != nbrPartitions.end(); it++ ) {
            for ( int idxParcel = 0; idxParcel < data->donePartitions[idxCluster].numParcels ; ++idxParcel ) {
                int parcelID = data->donePartitions[idxCluster].parcelIds[idxParcel];
                float disSample = this->GetDistanceHParcelAndPartition ( parcelID, *it );
                allDisNbr.push_back ( disSample );
            } // end for loop::idxParcel
        }

        vVarianceDis[idxCluster] = KML::StatsVariance ( allDisWithin ) /*KML::StatsMean(allDisWithin)*KML::StatsMean(allDisWithin)*/;
        allRatios[idxCluster] = KML::StatsMean ( allDisNbr ) /KML::StatsMean ( allDisWithin );
// 	if(400== data->donePartitions.size())
// 	{
// 	  cout<<allDisWithin<<endl;
// 	  cout<<allDisNbr<<endl;
// 	}
    } // end for loop::idxCluster




//   sigmaSquare=KML::StatsMean(vVarianceDis);
//   sigmaSquare=meanDis;
//   cout<<vMeanDis<<endl;
//   cout<<vVarianceDis<<endl;
//   cout<<KML::StatsMean(vMeanDis)<<" "<<KML::StatsMean(vVarianceDis)<<endl;

    //now estimate the likelihood ;
    //P(D|M,sigmaSquare) = multiplication of P(i | M, sigmaSquare);
    //P(i|M,sigmaSquare) = 1/N-nbrcluster*SUM(1/sqrt(2*pai*sigmaSquare)*exp(-1/2*Dis/sigmaSquare);

    double maxLikelihood = 0;
    for ( int idxCluster = 0; idxCluster < data->donePartitions.size() ; ++idxCluster ) {
        const std::set< int >& nbrPartitions = data->nbrPartitions[idxCluster];
        for ( int idxParcels = 0; idxParcels < data->donePartitions[idxCluster].numParcels ; ++idxParcels ) {

            int parcelID = data->donePartitions[idxCluster].parcelIds[idxParcels];
            double densityCurrentParcel = 0;

// 	    nbrs
            for ( std::set< int >::iterator it = nbrPartitions.begin(); it != nbrPartitions.end(); it++ ) {
                double sqrt2paisigmaS = sqrt ( 2*PAI*vVarianceDis[*it] );
                double invSqrt2paisigmaS = 1./sqrt2paisigmaS;
                double coef1= -0.5/vVarianceDis[*it];
                double coef2 = 1.f/nbrPartitions.size() *invSqrt2paisigmaS;
                float disSample = this->GetDistanceHParcelAndPartition ( parcelID, *it );
                densityCurrentParcel+= coef2*exp ( coef1*disSample*disSample );
            }

            //itself;
            {
                double sqrt2paisigmaS = sqrt ( 2*PAI*vVarianceDis[idxCluster] );
                double invSqrt2paisigmaS = 1./sqrt2paisigmaS;
                double coef1= -0.5/vVarianceDis[idxCluster];
                double coef2 = invSqrt2paisigmaS;
                float disSample = this->GetDistanceHParcelAndPartition ( parcelID, idxCluster );
                densityCurrentParcel+= coef2*exp ( coef1*disSample*disSample );
            }

            maxLikelihood+=log ( densityCurrentParcel );

        } // end for loop::idxParcels

    } // end for loop::idxCluster

//     maxLikelihood = log(maxLikelihood);


    aic = maxLikelihood- numberParameters;
    bic = maxLikelihood- 0.5* ( numberParameters ) *log ( data->numParcels );

    float part1ICL=0, part2ICL=0;
    for ( int idx = 0; idx < data->numParcels ; ++idx ) {
        part1ICL+= log ( idx+0.5*data->donePartitions.size() +1 );
    } // end for loop::idx

    for ( int idxCluster = 0; idxCluster < data->donePartitions.size() ; ++idxCluster ) {
        for ( int idxParcels = 0; idxParcels < data->donePartitions[idxCluster].parcelIds.size() ; ++idxParcels ) {
            part2ICL+= log ( idxParcels+1.5 );
        } // end for loop::idxParcels
    } // end for loop::idxCluster
    icl = bic - part1ICL + part2ICL;
    part1= maxLikelihood;
    cout<<"#clusters: "<<data->tobePartitions.size() <<" mean, sigma, part1,aic, bic, icl, ratio: "<<meanVar<<" "<<sigmaSquare<<" "<<part1<<" "<< aic<<" "<<bic<<" "<<icl<<" "<<KML::StatsMean ( allRatios ) <<endl;
}

//get the distance between the parcel and its partition;
float CFineParcellatorHCP::GetDistanceHParcelAndPartition ( int parcelIdx, int partitionIdx )
{
    //the distance dis=e^(-s/(d^p));
    //where s is the similarity;
    //d is the distance between the parcel and the partition;
    //p is the order for distance, the larger p is, the larger effect of distance is;

//   int parcelCenter = data->parcels[parcelIdx].centerVertexID;
//   int partitionCenter = data->donePartitions[partitionIdx].centerVertexID;
//   VectorType disEucVec = data->surf.GetPointCoords(parcelCenter) - data->surf.GetPointCoords(partitionCenter);
//   float disEuc = disEucVec.Norm();

//   cout<<data->parcels[parcelIdx].timeSeries<<endl;
//   cout<<data->donePartitions[partitionIdx].principalVector<<endl;
//     float sim = this->GetSimilarityTwoTimeSeries(data->parcels[parcelIdx].timeSeries, data->donePartitions[partitionIdx].principalVector);
//   cout<<data->parcels[parcelIdx].timeSeries<<endl;
//   cout<<data->donePartitions[partitionIdx].principalVector<<endl;
//     return 1-sim;
//     float result = atan2(0.1,sim-0.5);
//     return exp(-sim/*/disEuc*/);

    return this->GetDistanceH2TimeSeries ( data->parcels[parcelIdx].timeSeries, data->donePartitions[partitionIdx].principalVector );

} 

double CFineParcellatorHCP::GetDistanceH2TimeSeries ( std::vector< float >& s1, std::vector< float >& s2 )
{
    float sim = this->GetSimilarityTwoTimeSeries ( s1,s2 );
//     double result = atan2(0.1,sim-data->meanSimilarity);
    double result  =  1- sim;
    return result;
}

void CFineParcellatorHCP::DetermineMeanSimilarity ( int numSamples )
{
    int count=0;
    double allSims = 0.f;

    KML::CheckValueRange<> ( numSamples,100,data->numVertices );


    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc ( T );

    while ( count++ < numSamples ) {
        int parcelId1 = gsl_rng_uniform_int ( r,data->numParcels );
        int parcelId2 = gsl_rng_uniform_int ( r,data->numParcels );

        allSims+= this->GetSimilarityTwoTimeSeries ( data->parcels[parcelId1].timeSeries, data->parcels[parcelId2].timeSeries );
    }; //end of while;

    gsl_rng_free ( r );
    data->meanSimilarity = allSims/numSamples;
    cout<<"mean similarity of parcels with "<<numSamples<<" samples: "<<data->meanSimilarity<<endl;
}

void CFineParcellatorHCP::CalculateAreas4Partitions()
{
    data->surf.WarpPointLabel2CellLabel();
    data->allAreas.clear();
    data->allAreas.resize ( data->donePartitions.size() +1,0 );

    for ( int idxCell = 0; idxCell < data->surf.GetNumOfCells() ; ++idxCell ) {
        int cellLabel = data->surf.GetCellLabel ( idxCell );
        if ( -1!=cellLabel ) {
            data->allAreas[cellLabel]+= data->surf.GetCellArea ( idxCell );
        }
    } // end for loop::idxCell
    
    for ( int idx = 0; idx < data->donePartitions.size() ; ++idx )
    {
      data->donePartitions[idx].area = data->allAreas[idx]; 
    } // end for loop::idx

}


}//end of KML;






