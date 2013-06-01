/*
*	copyright by kaiming li (kaiming.li@emory.edu) 2012.
*	License : GPL v3
*
*/
#include <string>
#include <iostream>
#include "options.h"
#include "kaimingCommon.h"
#include <vector>
#include "kmlPartitions.h"
#include "armadillo"
#include "grangerCausality.h"
#include "tictoc.hpp"

using namespace Utilities;
using namespace KML;
using namespace std;
using namespace arma;




int main ( int argc, char **argv )
{

    string appName ( argv[0] );
    string appExample ( "lh.meta.CenterMetaInfo rh.meta.CenterMetaInfo labelMatch.lh.200.subjAtTmpl labelMatch.rh.200.subjAtTmpl outputbaseName[options]\
  \n\nDetails:\n\
  \tthis program calculates the connectivity pattern for the whole brain based on partitions.\n\
  \tlh.meta.CenterMetaInfo rh.meta.CenterMetaInfo contain all information for each partition;\n\
  \tIf need visualization on surface, additional arguments need to input; see usage for more information; \n\
  \tif granger causality used, then the postfix will be *.g; if pearson correlation is set, the postfix will be *.f. \n" );

    Option<bool> helpOpt ( string ( "-h" ),false,string ( "display this help information. " ),false,no_argument );
    Option<bool> grangerCausalityOpt ( "-g",false,"to calculate granger causality with postfix *.g, default: false", false,no_argument );
    Option<bool> pearsonOpt ( "-c",false, "to calculate the pearson correlation wich postfix *.f, default: false", false,no_argument );
    Option<string> surfaceOpt ( "-s","", "the surfaces lh and rh if visualization needed, no default", false, requires_2_arguments );
    Option<string> regSurfaceOpt ( "-r","", "the surfaces (registered to sphere) lh and rh if visualization needed, no default", false, requires_2_arguments );
    Option<float> pvalueOpt ( "-p",0.0001,"the pvalue to threshold  granger causality , default: 0.0001", false, requires_argument );
    Option<float> minCCOpt ( "-t",0.3, "the minimium pearson correlation, default: 0.3", false, requires_argument );
    Option<bool> verboseOpt ( "-v",false, "verboose mode if set", false, no_argument );
    Option<int> orderOpt ( "-o",5,"the order in the GC model, default 5", false,requires_argument );
    Option<bool> removeNbrOpt("-x",false,"remove the connectivity to neighbors, so that only long distance connectivity remains.",false,no_argument); 

    OptionParser cmdParser ( appName,appName+"\t"+appExample );
    cmdParser.add ( helpOpt );
    cmdParser.add ( grangerCausalityOpt );
    cmdParser.add ( pearsonOpt );
    cmdParser.add ( surfaceOpt );
    cmdParser.add ( regSurfaceOpt );
    cmdParser.add ( pvalueOpt );
    cmdParser.add ( minCCOpt );
    cmdParser.add ( verboseOpt );
    cmdParser.add ( orderOpt );
    cmdParser.add (removeNbrOpt); 

    if ( 6 > argc ) {
        cmdParser.usage();
        exit ( EXIT_FAILURE );
    }

    cmdParser.parse_command_line ( argc,argv,5 );
    if ( !grangerCausalityOpt.set() && !pearsonOpt.set() ) {
        cmdParser.usage();
        exit ( EXIT_FAILURE );
    }

    //the following 4 lines are for computational time of program;
    //put it where you want to start timing;
#ifdef TIC_TOC_HPP
    tictoc timer;
    timer.tic();
#endif




    vector<CPartition>  partitions,partitionsRh;
    fstream inStrm;
    KML::OpenReadStrmBinary ( inStrm,argv[1] ); //left parts;

    CPartition tmpPartition;
    while ( inStrm>>tmpPartition ) {
        partitions.push_back ( tmpPartition );
    }

    int lhSize = partitions.size();
    cout<<"lh parts size: "<<partitions.size() <<endl;


    inStrm.close();

    KML::OpenReadStrmBinary ( inStrm,argv[2] ); // right parts;
    while ( inStrm>>tmpPartition ) {
        partitions.push_back ( tmpPartition );
    }
    int rhSize = partitions.size()- lhSize;
    cout<<"rh parts size: "<<rhSize <<endl;
    inStrm.close();
    cout<<"number of vertices for first partition: "<< partitions[0].parcelIds.size()<<endl;

    int nAllParts = partitions.size();
    vector<int> allOrder;
    vector<int> aics,bics;

    vector<int> tmplMatchAtSubj; // get the subject counterpart of tmplate parts;

    {/**** Start of section::find matches; ****/
        vector<vector<int> > tmpLabelAtSubjlh, tmpLabelAtSubjrh;
        KML::ReadIntMatrix ( argv[3],tmpLabelAtSubjlh );
        KML::ReadIntMatrix ( argv[4],tmpLabelAtSubjrh );
        if ( tmpLabelAtSubjrh.size() != tmpLabelAtSubjlh.size() ) {
            cerr<<"Error:,unquall lh and rh sizes"<<endl;
            exit ( 1 );
        }
        for ( int idx = 0; idx < tmpLabelAtSubjlh.size() ; ++idx ) {
            tmplMatchAtSubj.push_back ( tmpLabelAtSubjlh[idx][0] ); // the first on of all possible matches;
        } // end for loop::idx

        for ( int idx = 0; idx < tmpLabelAtSubjrh.size() ; ++idx ) {
            tmplMatchAtSubj.push_back ( tmpLabelAtSubjrh[idx][0] +tmpLabelAtSubjlh.size() ); //+ offset from lh partitions; 
        } // end for loop::idx
    }/***** End  of section::find matches; ****/

//     cout<<tmplMatchAtSubj<<endl;


    fmat gcMat ( nAllParts,nAllParts );
    fmat fcMat ( nAllParts,nAllParts );

    if ( grangerCausalityOpt.set() ) {
        gcMat.fill ( 1 );
        int numSamples = partitions[0].principalVector.size();

        for ( int idx = 0; idx < nAllParts ; idx++ ) {
            if ( verboseOpt.set() && idx % 100 ==0 ) {
                cout<<"Calculate GC: "<< idx<<endl;
            }

            for ( int idx2 = idx+1; idx2 < nAllParts ; idx2++ ) {

                fmat data ( numSamples,2 );
                fvec x ( & partitions[idx].principalVector[0], numSamples,false );
                fvec y ( & partitions[idx2].principalVector[0], numSamples,false );
                data.col ( 0 ) = x;
                data.col ( 1 ) = y;
                CKmlGCA gca ( data, orderOpt.value() ); //using large order;

                const fmat& prob = gca.GCA();
                gcMat ( idx,idx2 ) = prob ( 0,1 );
                gcMat ( idx2,idx ) = prob ( 1,0 );
            } // end for loop::idx2


        } // end for loop::idx
        
         

        string outBaseName ( argv[5] );
        string outName= outBaseName;
        outName+=".g";
        gcMat.save ( outName, arma_binary );
        
 

        fmat threshGCMat = gcMat;
        for ( int x = 0; x < gcMat.n_rows ; ++x ) {
            for ( int y = 0; y < gcMat.n_cols ; ++y ) {
                if ( threshGCMat ( x,y ) < pvalueOpt.value() ) {
                    threshGCMat ( x,y ) = 1;
                } else {
                    threshGCMat ( x,y ) = 0;
                }
            } // end for loop::y
        } // end for loop::x
        outName+="t";
        threshGCMat.save ( outName, arma_binary );
        
        fmat gcFeatures = gcMat; 
        for ( int idx = 0; idx < tmplMatchAtSubj.size() ; ++idx )
        {
            int sidx = tmplMatchAtSubj[idx];
            for ( int idx2 = 0; idx2 < tmplMatchAtSubj.size() ; ++idx2 )
            {
                int sidx2 = tmplMatchAtSubj[idx2];
                gcFeatures(idx,idx2) = gcMat(sidx,sidx2); 
            } // end for loop::idx2
        } // end for loop::idx
        string featureName(argv[5]); 
        featureName+=".g.atTmplSpace"; 
        gcFeatures.save ( featureName, arma_binary );

    }

    if ( pearsonOpt.set() ) {
        fcMat.fill ( 0 );
        int numSamples = partitions[0].principalVector.size();

        for ( int idx = 0; idx < nAllParts ; idx++ ) {
            if ( verboseOpt.set() && idx % 100 ==0 ) {
                cout<<"Calculate FC: "<< idx<<endl;
            }
            for ( int idx2 = idx+1; idx2 < nAllParts ; idx2++ ) {
                float correl = KML::CoorelOfTwoSeries<> ( partitions[idx].principalVector,partitions[idx2].principalVector );
                fcMat ( idx,idx2 ) = correl;
                fcMat ( idx2,idx ) = correl;
            } // end for loop::idx2
        } // end for loop::idx
        
        //set diagnal to be max of all ; 
        for ( int idx = 0; idx < nAllParts ; ++idx )
	{
	  fcMat(idx,idx) = max(fcMat.col(idx)); 
	} // end for loop::idx
        
        if(removeNbrOpt.set())
	{
	  for ( int idxParts = 0; idxParts < partitions.size() ; ++idxParts )
	  {
	    int offset = idxParts< lhSize ? 0 : lhSize; 
	    vector<int>& nbrParts = partitions[idxParts].nbrPartitions;
	    for ( int idxNbr = 0; idxNbr < nbrParts.size() ; ++idxNbr )
	    {
	      int idxParts2 = nbrParts[idxNbr]+offset; 
	      fcMat(idxParts,idxParts2) = 0;
	      fcMat(idxParts2,idxParts) = 0;
	    } // end for loop::idxNbr
	    
	  } // end for loop::idxParts  
	} 
	
        string outBaseName ( argv[5] );
        string outName= outBaseName;
        outName+=".f";
        fcMat.save ( outName, arma_binary );

        fmat threshFCMat = fcMat;
        for ( int x = 0; x < fcMat.n_rows ; ++x ) {
            for ( int y = 0; y < fcMat.n_cols ; ++y ) {
                if ( abs ( threshFCMat ( x,y ) ) < minCCOpt.value() ) {
                    threshFCMat ( x,y ) = 0;
                }
            } // end for loop::y
        } // end for loop::x
        outName+="t";
        threshFCMat.save ( outName,arma_binary );
        
        
        fmat fcFeatures = fcMat; 
        for ( int idx = 0; idx < tmplMatchAtSubj.size() ; ++idx )
        {
            int sidx = tmplMatchAtSubj[idx];
            for ( int idx2 = 0; idx2 < tmplMatchAtSubj.size() ; ++idx2 )
            {
                int sidx2 = tmplMatchAtSubj[idx2];
                fcFeatures(idx,idx2) = fcMat(sidx,sidx2); 
            } // end for loop::idx2
        } // end for loop::idx
        string featureName(argv[5]); 
        featureName+=".f.atTmplSpace"; 
        fcFeatures.save ( featureName, arma_binary);

    }


    //the following 4 lines are for computational time of program;
    //put it where you want to end timing;
#ifdef TIC_TOC_HPP
    timer.toc();
    cout<<"CPU Time (ms) for program : "<<timer.totalTimeSec() <<endl;
#endif

     return 0;
} //end of main function;


