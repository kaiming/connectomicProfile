/*
*	copyright by kaiming li (kaiming.li@emory.edu) 2012.
*	License : GPL v3
*
*/
#include <string>
#include <iostream>
#include "options.h"
#include "kaimingCommon.h" 
#include "armadillo"
#include "apClustering.h"


using namespace Utilities;
using namespace KML;
using namespace std;
using namespace arma; 
 
 
int main(int argc, char **argv) {

  string appName(argv[0]); 
  string appExample("gRSNs commonRSNs idMaps \
	  	  	\n\nDetails:\n\
	                \tthis program finds the common RSNs from different group RSNs profiles\n\
	                \tcommonRSNs and idMaps are outputs for common profiles and RSN ID maps, respectively.\n\
	                \tgRSNsList is a file containing group RSN prfiles from different datasets. each line is a profile file \n\
		        \te.g.: \n\
		        \t\tfindCommonRsns comRsnPrfs idMaps gRSNsList\n"); 
  
  Option<bool> helpOpt(string("-h"),false,string("display this help information. "),false,no_argument);
  

  OptionParser cmdParser(appName,appName+"\t"+appExample);
  cmdParser.add(helpOpt);

  
  if(4 > argc){
    cmdParser.usage();
    exit(EXIT_FAILURE);
    
  }
 
  cmdParser.parse_command_line(argc,argv,3);
  vector<string> allgRSNs; 
  KML::ReadNameList(argv[1],allgRSNs); 
  
  
  vector<fmat> allGroupResults; 
  for ( int idxGroups = 0; idxGroups < allgRSNs.size() ; ++idxGroups )
  {
    fmat currentRst; 
    currentRst.load(allgRSNs[idxGroups]); 
    allGroupResults.push_back(currentRst);        
  } // end for loop::idxGroups
  
  
  if( 3 !=allGroupResults.size())
  {
    cout<<"not supported yet!"<<endl; 
    exit(1); 
  }
  
  
  //keep record of the id maps; 
  umat idRecs(allGroupResults[0].n_cols+allGroupResults[1].n_cols+ allGroupResults[2].n_cols, 2); //first col for group id, second col for prf id; 
  int idCount=0; 
  for ( int idxGroups = 0; idxGroups < allGroupResults.size() ; ++idxGroups )
  {
    for ( int idxPrf = 0; idxPrf < allGroupResults[idxGroups].n_cols ; ++idxPrf )
    {
      idRecs(idCount,0) = idxGroups;
      idRecs(idCount,1) = idxPrf; 
      idCount++; 
    } // end for loop::idxPrf
  } // end for loop::idxGroups
   
  
  list<int> v1, v2; 
  list<float> sims; 
  vector<float> allSims; 
  
  for ( int idx1 = 0; idx1 < idRecs.n_rows ; ++idx1 )
  {
    for ( int idx2 = idx1; idx2 < idRecs.n_rows ; ++idx2 )
    {
      	int g1 = idRecs(idx1,0); 
	int g2 = idRecs(idx2,0); 
	int r1 = idRecs(idx1,1); 
	int r2 = idRecs(idx2,1); 
      
      if(idRecs(idx1,0) !=idRecs(idx2,0)) //not in the same group; 
      {
	v1.push_back(idx1);
	v2.push_back(idx2);
	fmat correl = arma::cor( allGroupResults[g1].col(r1), allGroupResults[g2].col(r2) ); 
	sims.push_back(correl(0,0));	
      }
      
      fmat allCorr = arma::cor( allGroupResults[g1].col(r1), allGroupResults[g2].col(r2) ); 
      allSims.push_back(allCorr(0,0));
    } // end for loop::idx2
  } // end for loop::idx1
  
  CAPClustering apc; 
  apc.SetNPoints(idRecs.n_rows);
  apc.SetAplibName("/home/kli/bin/apclusterunix64.so");
  apc.Clustering(v1,v2,sims); 
  cout<<KML::StatsMean(allSims)<<endl;
  
  const  vector< std::vector< int > >& allClusters = apc.GetClusters(); 
  for ( int idxClusters = 0; idxClusters < allClusters.size() ; ++idxClusters )
  {

    vector<int> hasMem(3,0); 
    
    for ( int idxMem = 0; idxMem < allClusters[idxClusters].size() ; ++idxMem )
    {
      int id= allClusters[idxClusters][idxMem]; 
      int gid = idRecs(id,0); 
      int pid = idRecs(id,1); 
      hasMem[gid]++; 
      cout<<id<<" "; 
    } // end for loop::idxMem
    cout<<endl;
    if(hasMem[0]&&hasMem[1]&&hasMem[2])
          cout<<"common cluster, id: "<<idxClusters<<endl; 
  } // end for loop::idxClusters
  
 
 
  return 0;
} //end of main function;


