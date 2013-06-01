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
#include <indexer.h>
#include "cassert"
#include "omp.h"


using namespace Utilities;
using namespace KML;
using namespace std;
using namespace arma;

int main ( int argc, char **argv )
{

        string appName ( argv[0] );
        string appExample ( "templateSurf subjectSurf outputVids\
  \n\nDetails:\n\
  \tthis program finds for each vertex in templateSurf a matching vertex in subjectSurf that has minimum distance\n\
  \tthe input surfs should be spheres either in fs format or in vtk format;\n\
  \tthe outputVids is a binary array that has same length with templateSurf, with each element an ID in subjectSurf; \n\
  \te.g., appName  templ.sphere.reg lh.sphere.reg templ2sub.vid \n" );
        Option<bool> helpOpt ( string ( "-h" ),false,string ( "display this help information. " ),false,no_argument );
        Option<bool> isVTKOpt ( "-v",false, "whether or not in vtk format. default: false, i.e., in fs format; if set, then in vtk format", false,no_argument );
        Option<string> checkResultOpt ( "-r","check","the base name of output check surface: default check", false,requires_argument );


        OptionParser cmdParser ( appName,appName+"\t"+appExample );
        cmdParser.add ( helpOpt );
        cmdParser.add ( isVTKOpt );
        cmdParser.add ( checkResultOpt );


        if ( 4 > argc ) {
                cmdParser.usage();
                exit ( EXIT_FAILURE );
        }

        cmdParser.parse_command_line ( argc,argv,3 );


        CTriSurface surfTmpl ( argv[1],!isVTKOpt.set() );
        CTriSurface surfSubj ( argv[2],!isVTKOpt.set() );


        CIndexer indexer;

        /***
         * Start of section::initialize indexer;
         */
        {
                VectorType offset;
                offset.x= -200.;
                offset.y= -200.;
                offset.z= -200.;

                VectorType dims;
                dims.x= 4.;
                dims.y= 4.;
                dims.z= 4.;

                Vector3D<size_t> size ( 100,100,100 ); //default size is 1;

                indexer.SetDims ( dims );
                indexer.SetOffset ( offset );
                indexer.SetSize ( size );


        }
        /***
         * End of section::initialize indexer;
         */

        map<Vector3D<size_t>,vector<int> > mapIndexer2Vertices;
        /***
         * Start of section::build indexer for subjects;
         */
        {
                for ( int idxVSubj = 0; idxVSubj < surfSubj.GetNumOfPoints() ; ++idxVSubj ) {
                        Vector3D<size_t> grid = indexer.GetIndex ( surfSubj.GetPointCoords ( idxVSubj ) );
                        for ( int z = -1; z < 2 ; ++z ) {
                                for ( int y = -1; y < 2 ; ++y ) {
                                        for ( int x = -1; x < 2 ; ++x ) {
                                                Vector3D<size_t> newGrid ( grid.x+x, grid.y+y,grid.z+z );
                                                mapIndexer2Vertices[newGrid].push_back ( idxVSubj );
                                        } // end for loop::x

                                } // end for loop::y

                        } // end for loop::z

                } // end for loop::idxVSubj

        }
        /***
         * End of section::build indexer for subjects;
         */



        arma::uvec matches ( surfTmpl.GetNumOfPoints() );
        {/**** Start of section::find matches ****/
                #pragma omp parallel for
                for ( int idxVTmpl = 0; idxVTmpl < surfTmpl.GetNumOfPoints() ; ++idxVTmpl ) {

                        Vector3D<size_t> grid = indexer.GetIndex ( surfTmpl.GetPointCoords ( idxVTmpl ) );
                        vector<int>& candis = mapIndexer2Vertices[grid];
                        if ( candis.size() == 0 ) {
                                cout<<"ERROR: candis.size() = 0 "<<endl;
                                exit ( 1 );
                        }
                        arma::fvec allDis ( candis.size() );
                        allDis.fill ( 0 );
                        const VectorType& coordTarget = surfTmpl.GetPointCoords ( idxVTmpl );
                        for ( int idxCandi = 0; idxCandi < candis.size() ; ++idxCandi ) {
                                const VectorType& coorCandi = surfSubj.GetPointCoords ( candis[idxCandi] );
                                VectorType tmpDis = coorCandi - coordTarget;
                                allDis ( idxCandi ) = tmpDis.Norm();
                        } // end for loop::idxCandi

                        {/**** Start of section::find minimum distance ****/
                                uword minPoz ;
                                float minDis = allDis.min ( minPoz );
                                matches ( idxVTmpl ) = candis[minPoz];
                        }/***** End  of section::find minimum distance ****/

                } // end for loop::idxVTmpl
        }/***** End  of section::find matches ****/




        {/**** Start of section::write out  ****/
                matches.save ( argv[3],arma_binary );
        }/***** End  of section::write out  ****/


        if ( checkResultOpt.set() ) {
                /**** Start of section::check results ****/
                string baseName= checkResultOpt.value();
                surfTmpl.SaveAs ( baseName+".before.vtk" );
                for ( int idxV = 0; idxV < surfTmpl.GetNumOfPoints() ; ++idxV ) {
                        surfTmpl.ChangePointCoords ( idxV,surfSubj.GetPointCoords ( matches[idxV] ) );
                } // end for loop::idxV
                surfTmpl.SaveAs ( baseName+".after.vtk" );
        }/***** End  of section::check results ****/
  
        return 0;
} //end of main function;




