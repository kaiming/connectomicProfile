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
    string appExample ( "lh.parts.tmpl.vtk subj.lh.matchAt.tmpl tmpl.lh.matchAt.subj lh.parts.sub.vtk possibleMatches\
  \n\nDetails:\n\
  \tthis program finds the corresponding candidates of a template parcellation. \n\
  \tlh.parts.tmpl.vtk and lh.parts.sub.vtk should have same number of parts;\n\
  \tsubj.matchAt.tmpl is the match between subject vertex and tmpl vertex, nVertexInSub lines, each line is a vertex id in template; \n\
  \tpossibleMatches is a file with each line for possbile matches in subjects, with decreasing overlap with template " );


    Option<bool> helpOpt ( string ( "-h" ),false,string ( "display this help information. " ),false,no_argument );
    OptionParser cmdParser ( appName,appName+"\t"+appExample );    
    cmdParser.add ( helpOpt );

    if ( 6 > argc ) {
        cmdParser.usage();
        exit ( EXIT_FAILURE );
    }

    cmdParser.parse_command_line ( argc,argv,5 );

    CTriSurface tmplSurf ( argv[1] );
    CTriSurface subjSurf ( argv[4] );
    arma::uvec subjAtTmpl,tmplAtSubj;
    subjAtTmpl.load ( argv[2] );
    tmplAtSubj.load ( argv[3] );


    vector<int>& tmplLables = tmplSurf.GetAllPointLabel();
    vector<RGBTYPE>& tmplColors = tmplSurf.GetAllPointColor();
    vector<int>& subjLables = subjSurf.GetAllPointLabel();
    vector<RGBTYPE>& subjColors = subjSurf.GetAllPointColor();

    map<int, vector<int> > mapLabelVertices4Tmpl;

    {/**** Start of section::build label to vertices maps for template surface;  ****/
        for ( int idxVTmpl = 0; idxVTmpl < tmplSurf.GetNumOfPoints() ; ++idxVTmpl ) {
            int label = tmplLables[idxVTmpl];
            mapLabelVertices4Tmpl[label].push_back ( idxVTmpl );
        } // end for loop::idxVTmpl
    }/***** End  of section::build label to vertices maps for template surface;  ****/



    int numberOfLabelsTmpl = mapLabelVertices4Tmpl.size();
    vector<vector<int> > allMatches; 
  

    {/**** Start of section::main code ****/
        std::map< int, std::vector< int > >::iterator itmapLable = mapLabelVertices4Tmpl.begin();
        while ( itmapLable!=mapLabelVertices4Tmpl.end() ) {
            int lableT = itmapLable->first;
            if ( -1!=lableT ) {
                arma::ivec allPossibleLabels ( itmapLable->second.size() );
                allPossibleLabels.fill ( -2 ); //imporssible value;
                for ( int idxVCurrentLable = 0; idxVCurrentLable < itmapLable->second.size() ; ++idxVCurrentLable ) {
                    int correspondingVertexId = tmplAtSubj ( itmapLable->second[idxVCurrentLable] );
                    allPossibleLabels ( idxVCurrentLable ) = subjLables[correspondingVertexId];
                } // end for loop::idxVCurrentLable
                arma::uvec hc  = histc (allPossibleLabels, linspace<ivec>(-1,numberOfLabelsTmpl+10, numberOfLabelsTmpl+11 ));
                hc(0)=0; // remove medial wall; 
                
                arma::uword maxCount=0, maxIndex=-1; 
                maxCount = hc.max(maxIndex); 
                vector<int> currentPossibleMatches; 

                while(maxCount)
                {                   
                    currentPossibleMatches.push_back(-1+maxIndex);                  
                    hc(maxIndex)=0;
                    maxCount = hc.max(maxIndex);                     
                }
//                 for ( int idxCount = 1; idxCount < hc.n_elem ; ++idxCount )
//                 {
//                     if(hc(idxCount))
//                         currentPossibleMatches.push_back(-1+idxCount);
//                         
//                 } // end for loop::idxCount 

                allMatches.push_back(currentPossibleMatches);  
            }

            itmapLable++;
            
        }//end of while;
    }/***** End  of section::main code ****/
    
//     cout<<allMatches<<endl; 
    KML::SaveIntMatrix(argv[5],allMatches);
    
//     check 
//     map tmpl colors to subj; 
//     for ( int idxVSub = 0; idxVSub < subjSurf.GetNumOfPoints() ; ++idxVSub )
//     {
//         if(1==tmplLables[subjAtTmpl[idxVSub]])
//         {
//             subjSurf.SetPointColor(idxVSub,tmplColors[subjAtTmpl[idxVSub]]);
//         }
//         else
//             subjSurf.SetPointColor(idxVSub,RGBTYPE(0.2,0.2,0.2));
//         
//     } // end for loop::idxVSub
//     
//     subjSurf.SaveAs("check.vtk");
//     
//     for ( int idxVSub = 0; idxVSub < subjSurf.GetNumOfPoints() ; ++idxVSub )
//     {
//  
//         if(subjLables[idxVSub] == 11 )
//         {
//             cout<<idxVSub<<endl;
//             subjSurf.SetPointColor(idxVSub,RGBTYPE(1,0,0));
//         }
//         else if(subjLables[idxVSub] == 12 )
//         {
//             cout<<idxVSub<<endl;
//             subjSurf.SetPointColor(idxVSub,RGBTYPE(1,1,0));
//         }
//         else if (subjLables[idxVSub] == 25 )
//         {
//             cout<<idxVSub<<endl;
//             subjSurf.SetPointColor(idxVSub,RGBTYPE(1,0,1));
//         }
//         else if  (subjLables[idxVSub] == 26 )
//         {
//             cout<<idxVSub<<endl;
//             subjSurf.SetPointColor(idxVSub,RGBTYPE(0,1,0));
//         }
//         else
//             subjSurf.SetPointColor(idxVSub,RGBTYPE(0.2,0.2,0.2 ));
//         
//     } // end for loop::idxVSub
//     
//     subjSurf.SaveAs("check2.vtk");

    return 0;
} //end of main function;


