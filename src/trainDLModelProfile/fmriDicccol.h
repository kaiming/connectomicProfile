/*
*	copyright by kaiming li (kaiming.li@emory.edu) 2012.
*	License : GPL v3
*
*/
#ifndef __getfMRIDicccols__H
#define __getfMRIDicccols__H

#include <string> 
#include "armadillo"
#include "spams_dicts.h"
#include "spams_linalg.h"
#include "spams_decomp.h"

using SPAMS::Trainer; 
using SPAMS::ParamDictLearn; 
using SPAMS::Data; 
using SPAMS::SpMatrix;
using arma::fcube; 
using arma::fmat; 

namespace KML
{
  class CPrivacyfMRIDicccol; 
  class CfMRIDicccol
  {
  private: 

    CPrivacyfMRIDicccol* data; 
    void TrainingDL(const Data<float>& theData);
    void DecompLasso(const SPAMS::Matrix<float>& X,const SPAMS::Matrix<float>& D,SPAMS::SpMatrix<float>& spalpha);
    void SetTrainer(std::string DFile, const int batchsize = 256, const int NUM_THREADS=-1);
    void SetTrainer(std::string AFile, std::string BFile, std::string Dfile, const int itercount, const int batchsize, const int NUM_THREADS);



  public:
    
    CfMRIDicccol(); 
    ~CfMRIDicccol();
    
    void SetParamDicLearn(ParamDictLearn<float> param);
    void SetTrainer(const int k, const int batchsize = 256,const int NUM_THREADS=-1); 
    
    void SetData(std::string fileSeries); //read fmat series; 
    void SetData(fcube& allData); //read fcube data; 
    void SaveData(std::string fileName);  //save to fcube format; 

    const fmat& GetProfiles(void); 
    const fmat& GetProfileRatios(void); 
    
    int GetNROIs(void)const; 
    int GetNDims(void)const; 
    
    
    //set the file that contains directory of all subjects; e.g., ~/data2/fcp/zang/sub1/allReg; 
    void SetSubjectAllRegDirLists(string fileName); 
    //set the names of label matches, e.g., labelMatch.lh.resx.200.tmplAtSubj labelMatch.rh.resx.200.tmplAtSubj
    void SetNames4TmplLabelMatchAtSubj(string labelMatchlabelTmplAtSubj_LH, string labelMatchlabelTmplAtSubj_RH); 
    //set the names of partition centers at reg space, e.g., lh.parts.200.CentersonSphere.vtk rh.parts.200.CentersonSphere.vtk
    void SetNames4PartitionCenterSurfAtRegSpace(string Name4PartitionCenterSurfAtRegSpace_LH, string Name4PartitionCenterSurfAtRegSpace_RH); 
    //set the connectivity file name; 
    void SetConnectivityFileName(string connFileAtTmplSpace); 
    
    
    //save learned profile and coords models; 
    void SaveProfileAndCoordModels(string outputBaseName); 
    //save coordinate modesl; 
    void SaveCoordModels(string fileName); 
    //save trained profiles models; 
    void SaveProfileModels(std::string fileName); 
    
    
    //for clustering profiles; 
    float ComputeSims4Profiles(std::list< int >& v1, std::list< int >& v2, std::list< float >& sims, string simsFile = ""); 
    
    Trainer<float>* GetTrainer(void);
    ParamDictLearn<float>* GetParams(void); 
    
    //learns the connectivity profile for an ROI; 
    void LearnProfileModel4ROI(int idxROI);
    
    //learns the coordinate model for an ROI; i.e., mean and stdev of distribution; 
    void LearnCoordModel4ROI(int idxROI); 
    
    //learns the connectivity profiles for all ROIs; 
    void LearnProfileModels(void);
    
    //learns the anatomical models for all ROIs; note : must call learn profilemodels first; 
    void LearnCoordModels(void); 
    
    //learns both models; 
    void LearnProfileAndCoordModels(void); 
    
      
    
    
  };



}; 
#endif 

