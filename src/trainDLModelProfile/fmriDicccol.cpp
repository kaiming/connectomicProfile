/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012.
 *	License : GPL v3
 *
 */
#include "fmriDicccol.h"
#include "spams_dicts.h"
#include "kaimingCommon.h"
#include <triSurface.h>

using namespace KML;
using namespace std;
using namespace arma;

namespace KML
{
class CPrivacyfMRIDicccol
{

        friend class CfMRIDicccol;
private:

        Trainer<float>* trainer;
        ParamDictLearn<float>* params;
        SPAMS::Matrix<float> *X;
        fcube* allData;
        fmat* allProfiles;
        int K;
        int NUM_THREADS;
        int batch_size;
        fmat allProfileRatios;
        int nDims;
        int nROIs;
        int nSamples;
        vector<vector<int> > subjUsed4EachROI;
        vector<string> subjDirs;

        string lhLabelMatchName, rhLabelMatchName;
        string lhPartCenterAtRegName, rhPartCenterAtRegName;
        string connectivityFileName;
        vector<vector<VectorType> > allSampleCenters;
        vector<VectorType> coordMeanMdl4AllROIs;
        vector<float> coordDisMeanMdl4AllROIs;
        vector<float> coordDisStdevMdl4AllROIs;

public:
        CPrivacyfMRIDicccol() :
                trainer ( NULL ), params ( NULL ), X ( NULL ), allData ( NULL ), allProfiles ( NULL ),  K ( 5 ), NUM_THREADS ( -1 ), batch_size ( 30 ), lhLabelMatchName ( "labelMatch.lh.resx.200.tmplAtSubj" ),rhLabelMatchName ( "labelMatch.rh.resx.200.tmplAtSubj" ), lhPartCenterAtRegName ( "lh.parts.200.CentersonSphere.vtk" ), rhPartCenterAtRegName ( "rh.parts.200.CentersonSphere.vtk" ),connectivityFileName ( "connectivity.resx.200.f.atTmplSpace" ) {
        }
        ;

        ~CPrivacyfMRIDicccol() {
                if ( trainer )
                        delete trainer;
                if ( params )
                        delete params;
//         if (X)
//             delete X; //no need to delete since it used other structures data;
                if ( allData )
                        delete allData;
                if ( allProfiles )
                        delete allProfiles;
        }
        ;

};

CfMRIDicccol::CfMRIDicccol()
{
        data = new CPrivacyfMRIDicccol;
}

CfMRIDicccol::~CfMRIDicccol()
{
        delete data;
}

SPAMS::ParamDictLearn<float>* CfMRIDicccol::GetParams ( void )
{
        return data->params;
}
SPAMS::Trainer<float>* CfMRIDicccol::GetTrainer ( void )
{
        return data->trainer;
}
void CfMRIDicccol::SetParamDicLearn ( SPAMS::ParamDictLearn<float> param )
{
        data->params = new ParamDictLearn<float>();
        * ( data->params ) = param;
}
void CfMRIDicccol::SetTrainer ( const int k, const int batchsize, const int NUM_THREADS )
{
        data->K = k;
        data->batch_size = batchsize;
        data->NUM_THREADS = NUM_THREADS;
        if ( data->trainer )
                delete data->trainer;
        data->trainer = new SPAMS::Trainer<float> ( k, batchsize, NUM_THREADS );
}
void CfMRIDicccol::SetTrainer ( std::string DFile, const int batchsize, const int NUM_THREADS )
{

}
void CfMRIDicccol::SetTrainer ( std::string AFile, std::string BFile, std::string Dfile, const int itercount, const int batchsize, const int NUM_THREADS )
{

}
void CfMRIDicccol::SetData ( fcube& allData )
{
        data->allData = new fcube ( allData.memptr(), allData.n_rows, allData.n_cols, allData.n_slices );
        data->allProfiles = new fmat ( data->allData->n_cols, data->allData->n_cols );
}
void CfMRIDicccol::SetData ( std::string allRegDirsListFile )
{
        vector<string> allNames;
        KML::ReadNameList ( allRegDirsListFile, allNames );
	data->subjDirs= allNames; 
	
	for ( int idxNames = 0; idxNames < allNames.size() ; ++idxNames )
	{
	  allNames[idxNames]+="/";
	  allNames[idxNames]+= data->connectivityFileName; 
	} // end for loop::idxNames

        fmat currentfMat;
        currentfMat.load ( allNames[0] );

        int nRows = currentfMat.n_rows, nCols = currentfMat.n_cols, nSlices = allNames.size();
        data->allData = new fcube ( nRows, nCols, nSlices );
        data->allData->fill ( 0 );
        data->allData->slice ( 0 ) = currentfMat;
        for ( int idx = 1; idx < allNames.size(); ++idx ) {
                currentfMat.load ( allNames[idx] );
                data->allData->slice ( idx ) = currentfMat;
        } // end for loop::idx

        data->allProfiles = new fmat ( data->allData->n_cols, data->allData->n_cols );

        data->nROIs = nCols;
        data->nDims = nRows;
        data->nSamples = nSlices;
        data->subjUsed4EachROI.resize ( data->nROIs,vector<int>() );
        data->allSampleCenters.resize ( data->nROIs );
        data->coordMeanMdl4AllROIs.resize ( data->nROIs );
        data->coordDisMeanMdl4AllROIs.resize ( data->nROIs );
        data->coordDisStdevMdl4AllROIs.resize ( data->nROIs );

}

void CfMRIDicccol::SaveData ( std::string fileName )
{
        if ( data->allData )
                data->allData->save ( fileName );
        else {
                cerr << "Error: no data to save! " << endl;
                exit ( 1 );
        }
}

void CfMRIDicccol::TrainingDL ( const SPAMS::Data<float>& theData )
{
        if ( data->params ) {
                data->trainer->train ( theData, *data->params );
        } else {
                cerr << "Error, training parameters not set. Please set parameters first!  " << endl;
                exit ( 1 );
        }
}

void CfMRIDicccol::LearnProfileModel4ROI ( int idxROI )
{
        if ( data->trainer ) {
                if ( data->allData ) {
                        //get data ready for the training;
                        //check idxROI range first;
                        KML::CheckValueRange<int> ( idxROI, 0, data->allData->n_cols );

                        //then for each roi, there will be data->allData->n_slices samples;
                        fmat currentData ( data->allData->n_cols, data->allData->n_slices );
                        for ( int nSample = 0; nSample < data->allData->n_slices; ++nSample ) {
                                currentData.col ( nSample ) = data->allData->slice ( nSample ).col ( idxROI );
                        }
                        float* ptr = reinterpret_cast<float*> ( currentData.memptr() );
                        data->X = new SPAMS::Matrix<float> ( ptr, data->allData->n_cols, data->allData->n_slices );
                        if ( 65==idxROI ) {
                                string fileName ( "data." );
                                fileName+=idxROI;
                                currentData.save ( fileName );
                        }
                        //now training;
                        this->TrainingDL ( *data->X );
                } else {
                        cerr << "No Data available, please set data first!" << endl;
                        exit ( 1 );
                }
        } else {
                cerr << "Error, trainer not set. Please set trainer first!  " << endl;
                exit ( 1 );
        }
}
;

void CfMRIDicccol::DecompLasso ( const SPAMS::Matrix<float>& X, const SPAMS::Matrix<float>& D, SPAMS::SpMatrix<float>& spalpha )
{
        SPAMS::lasso2<float> ( X, D, spalpha, data->K, data->params->lambda, data->params->lambda2, data->params->mode, data->params->posAlpha, data->NUM_THREADS, NULL, 4 * data->K );
}

void CfMRIDicccol::LearnProfileModels ( void )
{
        if ( !data->allProfileRatios.n_cols )
                data->allProfileRatios = arma::randn<fmat> ( data->K,data->allData->n_cols );

        for ( int idxROI = 0; idxROI < data->allData->n_cols; ++idxROI ) {
                if ( 0==idxROI%50 ) {
//             cout<<"******************************************************"<<endl;
                        cout<<"Learning Profile Model for ROI Index : "<<idxROI<<endl;
                        cout<<"******************************************************"<<endl;
                }

                //need to reset trainer;
                this->SetTrainer ( data->K,data->batch_size,data->NUM_THREADS );

                this->LearnProfileModel4ROI ( idxROI );
                SPAMS::SpMatrix<float> spAlpha;
                SPAMS::Matrix<float> D, alpha;
                data->trainer->getD ( D );
                this->DecompLasso ( *data->X,D,spAlpha );
                spAlpha.toFull ( alpha );

                fmat alphaFull ( alpha.X(),alpha.m(),alpha.n() );

                //get the max of alpha;
                fvec sumAlpha = arma::sum ( alphaFull,1 );
//         if(65==idxROI)
// 	{
//             cout<<sumAlpha.t()<<endl;
//         }
                sumAlpha/=sum ( sumAlpha );
                arma::uword maxIdx;
                sumAlpha.max ( maxIdx );
                data->allProfileRatios.col ( idxROI ) = sumAlpha;


                //get max pos of alphaFull for coords model;
                vector<int> usedSubjs;
                for ( int idxCol = 0; idxCol < alphaFull.n_cols ; ++idxCol ) {
                        uword maxPos;
                        fvec currentCol = alphaFull.col ( idxCol );
                        float maxValue=currentCol.max ( maxPos );
                        if ( maxIdx == maxPos )
                                usedSubjs.push_back ( idxCol );
                } // end for loop::idxCol
                data->subjUsed4EachROI[idxROI] = usedSubjs;



                for ( int dim = 0; dim < data->allData->n_cols; ++dim ) {
                        data->allProfiles->operator() ( dim,idxROI ) = D ( dim,maxIdx );
                } //end for::dim

//         if(65==idxROI)
// 	{
//             fvec tmp = data->allProfiles->col(idxROI);
//             string fileName("dic.");
//             tmp.save(fileName+idxROI);
//         }

                if ( data->X )
                        delete data->X;
        } //end for::idxROI

}


void CfMRIDicccol::LearnCoordModels ( void )
{

        //prepare matched label;
        vector<vector<vector<int> >  > allMatches ( data->nSamples );
        {

                for ( int sid = 0; sid < data->nSamples ; ++sid ) {
                        vector<vector<int> > allMatchesSub;   // note first half and second half should use different surface file;
                        string lname = data->subjDirs[sid]+data->lhLabelMatchName;
                        string rname = data->subjDirs[sid]+data->rhLabelMatchName;
                        vector<vector<int> > rhMatches;
                        KML::ReadIntMatrix ( lname,allMatchesSub );
                        KML::ReadIntMatrix ( rname,rhMatches );

                        for ( int idx = 0; idx < rhMatches.size() ; ++idx ) {
                                allMatchesSub.push_back ( rhMatches[idx] );
                        } // end for loop::idx

                        allMatches[sid] = allMatchesSub;

                } // end for loop::sid
        }



        //prepare all surfs;
        vector<CTriSurface> allLH, allRH;
        for ( int sid = 0; sid < data->nSamples ; ++sid ) {
                string lname = data->subjDirs[sid]+data->lhPartCenterAtRegName;
                string rname = data->subjDirs[sid]+data->rhPartCenterAtRegName;
                CTriSurface lh ( lname ),rh ( rname );
                allLH.push_back ( lh );
                allRH.push_back ( rh );
        } // end for loop::sid

        for ( int idxROI = 0; idxROI < data->nROIs ; ++idxROI ) {
                bool isLeft = idxROI < data->nROIs/2 ? true : false;
                const vector<int>& usedSubjects = data->subjUsed4EachROI[idxROI];
                vector<VectorType> allCoords;

                for ( int sid = 0; sid < usedSubjects.size() ; ++sid ) {
                        vector<int>& allCorrespondences = allMatches[usedSubjects[sid]][idxROI]; //id is subject to half surf;
                        int pid = allCorrespondences[0]; // the first one overlapps most;
                        VectorType coord = isLeft ? allLH[usedSubjects[sid]].GetPointCoords ( pid ) :  allRH[usedSubjects[sid]].GetPointCoords ( pid );
                        allCoords.push_back ( coord );
                } // end for loop::sid

                data->allSampleCenters[idxROI] = allCoords;
                this->LearnCoordModel4ROI ( idxROI );
// 	string roiNames("center.samples.");
// 	roiNames+=idxROI;
// 	roiNames+=".vtk";
// 	KML::VisualizePointsUsingVerticeByVTK(allCoords,roiNames);
        } // end for loop::idxROI
}


void CfMRIDicccol::SaveCoordModels ( string fileName )
{
        KML::SaveStdVectorAsColumn<> ( fileName+".disMean", data->coordDisMeanMdl4AllROIs );
        KML::SaveStdVectorAsColumn<> ( fileName+".disStdev", data->coordDisStdevMdl4AllROIs );
        KML::SaveStdVectorAsColumn<> ( fileName+".mean",data->coordMeanMdl4AllROIs );
}


void CfMRIDicccol::SaveProfileModels ( std::string fileName )
{
        data->allProfiles->save ( fileName,arma_ascii );
        data->allProfileRatios.save ( fileName+".ratio",arma_ascii );
}
const fmat& CfMRIDicccol::GetProfiles ( void )
{
        if ( data->allProfiles )
                return * ( data->allProfiles );
        else {
                cerr<<"ERROR, please call LearnAllProfiles first"<<endl;
                exit ( 1 );
        }
}

const fmat& CfMRIDicccol::GetProfileRatios ( void )
{
        return data->allProfileRatios;
}
float CfMRIDicccol::ComputeSims4Profiles ( std::list< int >& v1, std::list< int >& v2, std::list< float >& sims, string simsFile )
{
        fmat correl = arma::cor ( * ( data->allProfiles ) );
        correl = abs ( correl );
        for ( int idxR = 0; idxR < correl.n_rows ; ++idxR ) {
                for ( int idxC = 0; idxC < correl.n_cols ; ++idxC ) {
                        if ( idxC!=idxR ) {
                                v1.push_back ( idxR );
                                v2.push_back ( idxC );
                                sims.push_back ( ( correl ( idxR,idxC ) ) );
                        } else
                                correl ( idxR,idxC ) =0;
                } // end for loop::idxC
        } // end for loop::idxR
        /*
            fstream simStrm;
            KML::OpenWriteStrmAscii(simStrm,"sims.txt");
            std::list< int >::iterator itv1=v1.begin();
            std::list< int >::iterator itv2=v2.begin();
            std::list< float >::iterator itsims=sims.begin();
            while(itv1!= v1.end())
            {
              simStrm<<*itv1+1<<" "<<*itv2+1<<" "<<*itsims<<endl;
              itv1++;itv2++;itsims++;

            }
            simStrm.close();
            exit(0);
            */

        if ( simsFile.size() )
                correl.save ( simsFile );

        //return mean value for preference;
        return mean ( mean ( ( correl ) ) );
}

int CfMRIDicccol::GetNDims ( void ) const
{
        return data->nDims;
}
int CfMRIDicccol::GetNROIs ( void ) const
{
        return data->nROIs;
}

void CfMRIDicccol::LearnProfileAndCoordModels ( void )
{
        this->LearnProfileModels();
        this->LearnCoordModels();
}

void CfMRIDicccol::LearnCoordModel4ROI ( int idxROI )
{
        //get mean location;
        vector<VectorType>& currCenters = data->allSampleCenters[idxROI];
        //normalize the coordinates to same radius, and convert euclidean coords to spherical coords;
        float expectedRadius = 100.f;
        for ( int idx = 0; idx < currCenters.size() ; ++idx ) {
                float ratio = currCenters[idx].Norm() /expectedRadius;
                currCenters[idx]/=ratio;
                VectorType sphCorrd;
                KML::ConvertCoordsEuclidean2Spherical ( currCenters[idx],sphCorrd );
                currCenters[idx] = sphCorrd;
        } // end for loop::idx


        VectorType meanLoc;
        for ( int idx = 0; idx < currCenters.size() ; ++idx ) {
                meanLoc+=currCenters[idx];
        } // end for loop::idx
        meanLoc/=currCenters.size();
        data->coordMeanMdl4AllROIs[idxROI] = meanLoc;

        //get distance;
        fvec allDis ( currCenters.size() );
        allDis.fill ( 0 );
        for ( int idx = 0; idx < currCenters.size() ; ++idx ) {
                allDis[idx] = KML::GetGeodesicDisSphereicalCoords ( meanLoc,currCenters[idx] );
        } // end for loop::idx
        data->coordDisMeanMdl4AllROIs[idxROI] = mean ( allDis );
        data->coordDisStdevMdl4AllROIs[idxROI] = stddev<> ( allDis );

}

void CfMRIDicccol::SetSubjectAllRegDirLists ( string fileName )
{
        KML::ReadNameList ( fileName,data->subjDirs );
        if ( data->subjDirs.size() != data->nSamples ) {
                cerr<<" Error: unqual size of samples/subjects. "<<endl;
                exit ( 1 );
        }
}

void CfMRIDicccol::SetNames4PartitionCenterSurfAtRegSpace ( string Name4PartitionCenterSurfAtRegSpace_LH, string Name4PartitionCenterSurfAtRegSpace_RH )
{
        data->lhPartCenterAtRegName = Name4PartitionCenterSurfAtRegSpace_LH;
        data->rhPartCenterAtRegName = Name4PartitionCenterSurfAtRegSpace_RH;
}
void CfMRIDicccol::SetConnectivityFileName ( string connFileAtTmplSpace )
{
        data->connectivityFileName = connFileAtTmplSpace;

}

void CfMRIDicccol::SetNames4TmplLabelMatchAtSubj ( string labelMatchlabelTmplAtSubj_LH, string labelMatchlabelTmplAtSubj_RH )
{
        data->lhLabelMatchName = labelMatchlabelTmplAtSubj_LH;
        data->rhLabelMatchName = labelMatchlabelTmplAtSubj_RH;
}

void CfMRIDicccol::SaveProfileAndCoordModels ( string outputBaseName )
{
        this->SaveProfileModels ( outputBaseName+".prf" );
        this->SaveCoordModels ( outputBaseName+".coord" );

}


}//end of workspace KML:
;

