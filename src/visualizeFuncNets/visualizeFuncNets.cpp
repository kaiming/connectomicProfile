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
#include "armadillo"



using namespace Utilities;
using namespace KML;
using namespace std;
using namespace arma;


int main ( int argc, char **argv )
{

    string appName ( argv[0] );
    string appExample ( "tmpl.parts.lh.vtk tmpl.parts.rh.vtk clusterProfiles outputBase\
			\n\nDetails:\n\
			\tthis program generates files for visualization given profiles of RSNs\n\
			\tthe profiles must be expressed in the tmplate parcellation scheme, or it will mess up. ");

    Option<bool> helpOpt ( string ( "-h" ),false,string ( "display this help information. " ),false,no_argument );
    Option<int> idxOpt("-n",14,"the roi whose networks will be written out. if set to -1, then all will be written out.", false,requires_argument);
    Option<float> shiftTmplSurfOpt("-s",0,"the distance to seperate left and right semisphere, default: 0", false,requires_argument); 
    Option<float> profileThreshOpt("-t",0.001,"the threshold of profile , value under the threshold will not be considered, default: 0.001", false,requires_argument); 

    OptionParser cmdParser ( appName,appName+"\t"+appExample );
    cmdParser.add ( helpOpt );
    cmdParser.add( idxOpt );
    cmdParser.add(shiftTmplSurfOpt);
    cmdParser.add(profileThreshOpt);

    if ( 5 > argc ) {
        cmdParser.usage();
        exit ( EXIT_FAILURE );
    }

    cmdParser.parse_command_line ( argc,argv,4 );

    CTriSurface lh ( argv[1] ),rh ( argv[2] );
    
    if(shiftTmplSurfOpt.set())
    {
      for ( int idxV = 0; idxV < lh.GetNumOfPoints() ; ++idxV )
      {
	VectorType coord= lh.GetPointCoords(idxV);
	coord.x += shiftTmplSurfOpt.value(); 
	lh.ChangePointCoords(idxV, coord);
      } // end for loop::idxV      
    } 
    fmat clusterProfiles;
    clusterProfiles.load ( argv[3] );
 
 
    int netID=0;
    vector<vector<float> > allAtrri;

    //write out the vertices alone;
    vector<vector<int> > vidLh,vidRH;
    vector<vector<float> > dataLh,dataRh;
    
    int numPartitions = clusterProfiles.n_rows; 
    
    for ( int idxRSN = 0; idxRSN < clusterProfiles.n_cols ; ++idxRSN )
    {
        subview_col< float > currentDic= clusterProfiles.col(idxRSN); 
      	
        vector<float> attriLH ( lh.GetNumOfPoints(),0 ),attriRH ( rh.GetNumOfPoints(),0 );  
        vector<int> tmpVidLh,tmpVidRh;
        vector<float> tmpDataLh,tmpDataRh;
	
        for ( int partIdx = 0; partIdx < numPartitions ; ++partIdx ) {
            if ( partIdx<numPartitions/2 && abs ( currentDic ( partIdx ) ) > profileThreshOpt.value() ) {
                for ( int idxV = 0; idxV < lh.GetNumOfPoints() ; ++idxV ) {
                    if ( partIdx == lh.GetPointLabel ( idxV ) ) {
                        attriLH[idxV] = currentDic ( partIdx );
                        tmpVidLh.push_back ( idxV );
                        tmpDataLh.push_back ( currentDic ( partIdx ) );
                    }
                } // end for loop::idxV
            } else if ( partIdx >= numPartitions/2  && abs ( currentDic ( partIdx ) ) >profileThreshOpt.value() ) {
                for ( int idxRh = 0; idxRh < rh.GetNumOfPoints() ; ++idxRh ) {
                    if ( partIdx == rh.GetPointLabel ( idxRh ) + numPartitions/2 ) {
                        attriRH[idxRh] = currentDic ( partIdx );
                        tmpVidRh.push_back ( idxRh );
                        tmpDataRh.push_back ( currentDic ( partIdx ) );
                    }
                } // end for loop::idxRh
            }

        } // end for loop::partIdx

        string attriName ( "net." );
        attriName+=netID;
        lh.AddPointDataScalar ( attriLH,attriName );
        rh.AddPointDataScalar ( attriRH,attriName );



        vector<float> tmpAttri;
        tmpAttri.insert ( tmpAttri.end(),attriLH.begin(),attriLH.end() );
        tmpAttri.insert ( tmpAttri.end(),attriRH.begin(),attriRH.end() );
        allAtrri.push_back ( tmpAttri );

        vidLh.push_back ( tmpVidLh );
        vidRH.push_back ( tmpVidRh );
        dataLh.push_back ( tmpDataLh );
        dataRh.push_back ( tmpDataRh );

        netID++;
      
    } // end for loop::idxRSN
 
 
    cout<<"write all networks as attributes to : "<<string(argv[4])+".allRSNs.vtk"<<endl;
    lh.SaveAs(string(argv[4])+".lh.allRSNs.vtk");
    rh.SaveAs(string(argv[4])+".rh.allRSNs.vtk");
    lh.CombineSurf ( rh );
    lh.SaveAs ( string(argv[4])+".allRSNs.vtk");
 

    if(idxOpt.set())
    {
        cout<<"write individual functional networks ...";
        for ( int idxCluter = 0; idxCluter < vidLh.size() ; ++idxCluter ) {
            if (  -1==idxOpt.value() || idxCluter==idxOpt.value() ) {
                vector<VectorType> allVids;
                vector<float> allAttriData;
                for ( int idxL = 0; idxL < vidLh[idxCluter].size() ; ++idxL ) {
                    allVids.push_back ( lh.GetPointCoords ( vidLh[idxCluter][idxL] ) );
                    allAttriData.push_back ( dataLh[idxCluter][idxL] );
                } // end for loop::idxL
                string fileName ( string(argv[4])+".RSN." );
                fileName+=idxCluter;
                fileName+=".lh.vtk";
                KML::VisualizePointsUsingVerticeByVTKWithPointData ( allVids,allAttriData,fileName );
            }
            if ( -1==idxOpt.value() || idxCluter==idxOpt.value() ) {
                vector<VectorType> allVids;
                vector<float> allAttriData;
                for ( int idxR = 0; idxR < vidRH[idxCluter].size() ; ++idxR ) {
                    allVids.push_back ( rh.GetPointCoords ( vidRH[idxCluter][idxR] ) );
                    allAttriData.push_back ( dataRh[idxCluter][idxR] );
                } // end for loop::idxR
                string fileName ( string(argv[4])+".RSN." );
                fileName+=idxCluter;
                fileName+=".rh.vtk";
                KML::VisualizePointsUsingVerticeByVTKWithPointData ( allVids,allAttriData,fileName );
            };
            cout<<"done"<<endl;

        } // end for loop::idxCluter
    }//end of if;
 

    return 0;
} //end of main function;


//original core codes; 
//     CTriSurface lh ( argv[1] ),rh ( argv[2] );
//     ivec idx;
//     fmat sims;
//     fmat allDics;
// 
//     idx.load ( argv[3] );
//     allDics.load ( argv[4] );
//     sims.load ( argv[5] );
// 
//     if ( idx.n_elem != allDics.n_cols ||idx.n_elem != sims.n_cols ) {
//         cerr<<"Unequal size of idx/allDics/sims, ERROR! "<<endl;
//         exit ( 1 );
//     }
// 
//     vector<int> centers;
//     map<int, vector<int> > nets;
//     map<int, fvec> profiles;
// 
//     {   /**** Start of section::get examplars of clustering; ****/
//         for ( int index = 0; index < idx.n_elem ; ++index ) {
//             if ( index == ( int ) idx ( index ) ) {
//                 centers.push_back ( index );
//             }
//             nets[ ( int ) idx ( index ) ].push_back ( index );
//             if ( profiles[ ( int ) idx ( index )].n_elem==0 ) {
//                 profiles[ ( int ) idx ( index )]=allDics.col ( index );
//             } else {
//                 profiles[ ( int ) idx ( index )]+=allDics.col ( index );
//             }
//         } // end for loop::idx
//         cout<<"centers: " <<centers.size() <<endl<<centers<<endl;
//     }/***** End  of section::get examplars of clustering; ****/
// 
// 
// 
//     std::map< int, std::vector< int > >::iterator it= nets.begin();
//     int netID=0;
//     vector<vector<float> > allAtrri;
// 
//     //write out the vertices alone;
//     vector<vector<int> > vidLh,vidRH;
//     vector<vector<float> > dataLh,dataRh;
//     while ( it!=nets.end() ) {
// 
// 
//         int currExamplar = it->first;
//         vector<int>& members = it->second;
//         profiles[currExamplar]/=members.size(); //get average;
// 
//         vector<float> attriLH ( lh.GetNumOfPoints(),0 ),attriRH ( rh.GetNumOfPoints(),0 );
//         fvec& currentDic= profiles[currExamplar];
// 
// 
//         vector<int> tmpVidLh,tmpVidRh;
//         vector<float> tmpDataLh,tmpDataRh;
//         for ( int partIdx = 0; partIdx < idx.n_elem ; ++partIdx ) {
//             if ( partIdx<idx.n_elem/2 && abs ( currentDic ( partIdx ) ) >0.001 ) {
//                 for ( int idxV = 0; idxV < lh.GetNumOfPoints() ; ++idxV ) {
//                     if ( partIdx == lh.GetPointLabel ( idxV ) ) {
//                         attriLH[idxV] = currentDic ( partIdx );
//                         tmpVidLh.push_back ( idxV );
//                         tmpDataLh.push_back ( currentDic ( partIdx ) );
//                     }
//                 } // end for loop::idxV
//             } else if ( partIdx >= idx.n_elem/2  && abs ( currentDic ( partIdx ) ) >0.001 ) {
//                 for ( int idxRh = 0; idxRh < rh.GetNumOfPoints() ; ++idxRh ) {
//                     if ( partIdx == rh.GetPointLabel ( idxRh ) + idx.n_elem/2 ) {
//                         attriRH[idxRh] = currentDic ( partIdx );
//                         tmpVidRh.push_back ( idxRh );
//                         tmpDataRh.push_back ( currentDic ( partIdx ) );
//                     }
//                 } // end for loop::idxRh
//             }
// 
//         } // end for loop::partIdx
// 
//         string attriName ( "net." );
//         attriName+=netID;
//         lh.AddPointDataScalar ( attriLH,attriName );
//         rh.AddPointDataScalar ( attriRH,attriName );
// 
// 
// 
//         vector<float> tmpAttri;
//         tmpAttri.insert ( tmpAttri.end(),attriLH.begin(),attriLH.end() );
//         tmpAttri.insert ( tmpAttri.end(),attriRH.begin(),attriRH.end() );
//         allAtrri.push_back ( tmpAttri );
// 
//         vidLh.push_back ( tmpVidLh );
//         vidRH.push_back ( tmpVidRh );
//         dataLh.push_back ( tmpDataLh );
//         dataRh.push_back ( tmpDataRh );
// 
// 
//         it++;
//         netID++;
//     }
//     cout<<"write all networks as attributes to : funcNet.whole.vtk"<<endl;
//     lh.CombineSurf ( rh );
//     lh.SaveAs ( "funcNet.whole.vtk" );