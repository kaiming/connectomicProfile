/*
 *	copyright by kaiming li (kaiming.li@emory.edu) 2012.
 *	License : GPL v3
 *
 */
#include <string>
#include <iostream>
#include "options.h"
#include "kaimingCommon.h"
#include "kmlPartitions.h"
#include "armadillo"
#include "vector"
#include "triSurface.h"
#include "algorithm"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_rng.h"

using namespace Utilities;
using namespace KML;
using namespace std;
using namespace arma;



int main ( int argc, char **argv )
{

    string appName ( argv[0] );
    string appExample ( "mdlConfigFile parcelConfigFile subAllRegDir outBaseName \n\
		      \n\nDetails:\n\
		      \tthis program finds individual functional networks based on anatomical and profile model information\n\
		      \tmdlConfigFile contains the full-path  model files, including (in order): \n\
		      \t\t\tcluster.idx cluster.sims coord.disMean coord.disStdev coord.mean prf prf.cluster prf.ratio \n\
		      \tparcelConfigFile is a file specifying the names of relavant inputfiles, including (in order):  \n\
		      \t\t\tlh.CenterMetaInfo lh.vtk rh.CenterMetaInfo rh.vtk subj.lh.matchAt.tmpl subj.rh.matchAt.tmpl tmpl.lh.matchAt.subj tmpl.rh.matchAt.subj connectivity.atTmplSpace. \n\
		      \tthe energy function used : \n\
	              \tCi=Min(lamda*IntF(Ci,Ti)+(1-lamda)*ExtF(Ci,Ti)) \n\
	              \tNote, all outputs will be written to subAllRegDir" );

    Option<bool> helpOpt ( string ( "-h" ), false, string ( "display this help information. " ), false, no_argument );
    Option<float> lamdaOpt ( "-l", 0.1, "the weight of anatomical information, default: 0.1", false, requires_argument );
    Option<bool> writeIndivFuncNetOpt ( "-s", false, "whether or not write out the individual functional networks as vtk files, default: false", false, no_argument );
    Option<float> minThreshProfileOpt ( "-e", 0.01, "the minimal values in a profile to be considered, default: 0.01", false, requires_argument );

    OptionParser cmdParser ( appName, appName + "\t" + appExample );
    cmdParser.add ( helpOpt );
    cmdParser.add ( lamdaOpt );
    cmdParser.add ( writeIndivFuncNetOpt );
    cmdParser.add ( minThreshProfileOpt );

    if ( 5 > argc ) {
        cmdParser.usage();
        exit ( EXIT_FAILURE );
    }

    cmdParser.parse_command_line ( argc, argv, 4 );




    fmat centerPrfs;
    vector<VectorType> coordMdl;
    vector<float> coordDisStdev, coordDisMean;
    fmat connectivityPatternAtTmplSpace,connectivityPatternAtSubjSpace;
    vector<string> allMdlFiles;
    vector<vector<int> > roiLabels;
    KML::ReadNameList ( argv[1],allMdlFiles );
    {   /**** Start of section::read coords models; processing argv[1] ****/

        centerPrfs.load ( allMdlFiles[6] ); // read center profiles; the 7th in the list;

        fstream inStrm;
        KML::OpenReadStrmAscii ( inStrm, allMdlFiles[4] ); // the roi coords at template space ; the 5th;
        VectorType tmpCoord;
        while ( inStrm >> tmpCoord.x >> tmpCoord.y >> tmpCoord.z && tmpCoord.x != 0 ) {
            coordMdl.push_back ( tmpCoord );
        }

        KML::ReadColumnVector<> ( allMdlFiles[2], coordDisMean ); //dis mean; the 3rd
        KML::ReadColumnVector<> ( allMdlFiles[3], coordDisStdev ); //dis stdev; the 4th;
        KML::ReadIntMatrix(allMdlFiles[10],roiLabels);

    }/***** End  of section::read coords models; processing argv[1] ****/


    vector<vector<int> > labelMatchSubjAtTmpl, labelMatchTmplAtSubj; // for label matches, first lh, and then rh;
    vector<CPartition> subjParts; // for all partitions; first lh, and then rh;
    string subjectRootDir ( argv[3] );
    vector<string> parcelConfigFileName;
    KML::ReadNameList ( argv[2],parcelConfigFileName );
    subjectRootDir+="/";
    string logName(subjectRootDir+string(argv[4])+".log");
    fstream outStrm;
    KML::OpenWriteStrmAscii(outStrm,logName);

    {   /******** start of section::read parcel information; processing argv[2] ****/

        //labelMatchSubjAtTmpl
        vector<vector<int> > rightMatch;
        KML::ReadIntMatrix ( subjectRootDir+parcelConfigFileName[4], labelMatchSubjAtTmpl ); //the 5th one; subj.lh.matchAt.tmpl
        KML::ReadIntMatrix ( subjectRootDir+parcelConfigFileName[5], rightMatch ); // the sixth one;
        for ( int idx = 0; idx < rightMatch.size() ; ++idx ) {
            labelMatchSubjAtTmpl.push_back ( rightMatch[idx] );
        } // end for loop::idx

        //labelMatchTmplAtSubj
        rightMatch.clear();
        KML::ReadIntMatrix ( subjectRootDir+parcelConfigFileName[6], labelMatchTmplAtSubj ); //the 7th one; tmpl.lh.matchAt.subj
        KML::ReadIntMatrix ( subjectRootDir+parcelConfigFileName[7], rightMatch ); // the sixth one;
        for ( int idx = 0; idx < rightMatch.size() ; ++idx ) {
            labelMatchTmplAtSubj.push_back ( rightMatch[idx] );
        } // end for loop::idx



        {   /**** Start of section::read subject parts;  ****/

            fstream instrm;
            KML::OpenReadStrmBinary ( instrm, subjectRootDir+parcelConfigFileName[0] ); //subject parcellations::LH; // the first one in the config file;
            CPartition currentParts;
            while ( instrm >> currentParts ) {
                subjParts.push_back ( currentParts );
            }
            instrm.close();

            KML::OpenReadStrmBinary ( instrm, subjectRootDir+parcelConfigFileName[2] ); //subject parcellations::LH;
            while ( instrm >> currentParts ) {
                subjParts.push_back ( currentParts );
            }
            instrm.close();

        }/***** End  of section::read  subject parts;  ****/

        connectivityPatternAtTmplSpace.load ( subjectRootDir+parcelConfigFileName[8] ); // for connectivity pattern at tmpl space ;
        connectivityPatternAtSubjSpace.load ( subjectRootDir+parcelConfigFileName[9] ); // for connectivity pattern at subj space;

    }/******** end of section::read parcel information; processing argv[2] ****/


    fmat allExtEng= centerPrfs;
    allExtEng.fill ( 0 );
    fmat allIntEng = allExtEng;
    fmat allTotalEng = allExtEng;  //keep records of interal  external and total energy for all involved nodes of all func nets;
    fmat allCorrelation = allExtEng;

//     fmat allIndividualFuncNetProfiles = centerPrfs;     //keep record all the profiles at tmplate space;
    fmat allIndividualFuncNetsNodes ( centerPrfs.n_rows,centerPrfs.n_cols ); //keep record of all involved partitions;
    allIndividualFuncNetsNodes.fill ( -1.f );
//     allIndividualFuncNetProfiles.fill(0);
    fmat mapSub2tmpIdx = centerPrfs;  // for each net, map the subj region idx to tmp region idx;
    mapSub2tmpIdx.fill ( -1 );


    //main loops;
    int stepSize = 0;

    //cout << "TmplIdx | SubIdx | IntEnergy | ExtEnergy | TotalEnergy " << endl;
    //cout << fixed;
    outStrm<< "TmplIdx | SubIdx | IntEnergy | ExtEnergy | TotalEnergy " << endl;
    outStrm<< fixed;

    for ( int idxNets = 0; idxNets < centerPrfs.n_cols ; ++idxNets ) {
//         cout<<"############################# Processing network idx :"<< idxNets<<" ##############################################"<<endl;
        outStrm<<"############################# Processing network idx :"<< idxNets<<" ##############################################"<<endl;
//         cout<<roiLabels[idxNets]<<endl;



        for ( int idxROI = 1; idxROI < roiLabels[idxNets].size(); ++idxROI ) {


            int idxTmplRegions = roiLabels[idxNets][idxROI];
            if (1 ) { // will only consider significant regions;
                int step = idxTmplRegions < centerPrfs.n_rows/2 ? 0 : centerPrfs.n_rows/2;
                int tmplateLabel = idxTmplRegions - step;
                int subMatchedLabel = labelMatchTmplAtSubj[idxTmplRegions][0];
                int subPartIdx = step+ subMatchedLabel;

//                 //cout<<idxTmplRegions<<"/"<<subPartIdx<<"/"<<step<<"/"<<tmplateLabel<<"/"<<subMatchedLabel<<"/"<<centerPrfs(idxTmplRegions,idxNets)<<endl;

                vector<int> allCandiRegionsVec = subjParts[subPartIdx].nbrPartitions;
                allCandiRegionsVec.push_back ( subMatchedLabel );

                //internal energies are the same for all candidates, but external energy is network specific;
                fvec candiExtEng ( allCandiRegionsVec.size() ), candiTotalEng ( allCandiRegionsVec.size() ),candiIntEng ( allCandiRegionsVec.size() );
                candiExtEng.fill ( 1.e10 );
                candiTotalEng.fill ( 1.e10 );
                candiIntEng.fill ( 1.e10 );

                VectorType& mdlSphereCoord = coordMdl[idxTmplRegions];
                {   /**** Start of section::internal energy ****/
                    for ( int idxCandi = 0; idxCandi < allCandiRegionsVec.size(); ++idxCandi ) {
                        int labelCandi = allCandiRegionsVec[idxCandi];
                        VectorType& sphereCoord = subjParts[step+labelCandi].sphereCoords;
                        sphereCoord.x= mdlSphereCoord.x;
                        float geoDis = KML::GetGeodesicDisSphereicalCoords ( sphereCoord, mdlSphereCoord );
                        float value = ( geoDis - coordDisMean[idxTmplRegions] ) / ( 3 * coordDisStdev[idxTmplRegions] );
                        candiIntEng[idxCandi] = value;
                    } // end for loop::idxCandi
                }/***** End  of section::internal energy ****/

                {   /**** Start of section::ext energy ****/
                    subview_col<float> prfTi = centerPrfs.col ( idxNets );

                    uvec posValueIdx = find ( prfTi > minThreshProfileOpt.value() );
                    uvec negValueIdx = find ( prfTi < -minThreshProfileOpt.value() );

                    for ( int idxCandi = 0; idxCandi < allCandiRegionsVec.size(); ++idxCandi ) {
                        int candiLabelSubj = allCandiRegionsVec[idxCandi];
                        int tmplateCorrepondingIdx = labelMatchSubjAtTmpl[step+candiLabelSubj][0] + step;
                        subview_col<float> prfCi = connectivityPatternAtTmplSpace.col ( tmplateCorrepondingIdx );

                        fvec sortedPrfCi = arma::sort ( prfCi, 1 ); //decending sort;
                        float posThreshValue = sortedPrfCi[posValueIdx.n_elem];
                        float negThreshValue = sortedPrfCi[sortedPrfCi.n_elem - 1 - negValueIdx.n_elem];

                        for ( int idxPrfCi = 0; idxPrfCi < prfCi.n_elem; ++idxPrfCi ) {
                            prfCi[idxPrfCi] = ( prfCi[idxPrfCi] > posThreshValue || prfCi[idxPrfCi] < negThreshValue ) ? prfCi[idxPrfCi] : 0;
                        } // end for loop::idxPrfCi
                        fmat correls = arma::cor ( prfCi, prfTi );
                        candiExtEng[idxCandi] = 1 - sgn ( centerPrfs ( idxTmplRegions,idxNets ) ) * ( correls ( 0, 0 ) );
                    } // end for loop::idxCandi
                }/***** End  of section::ext energy ****/

                {   /**** Start of section::total engergy ****/
                    float lamda =0.f;
                    if(idxNets == 7 || idxNets == 13 || idxNets == 23 )
                        lamda = 1.f;
                    else
                        lamda = lamdaOpt.value();
                    candiTotalEng = lamda * candiIntEng + ( 1.f - lamda ) * candiExtEng;
                }/***** End  of section::total engergy ****/

                //find minimal engery now;
                uword minPox;
                float minValue;
                minValue = candiTotalEng.min ( minPox );

//                 cout<<"Candis for "<<roiLabels[idxNets][idxROI]<<" : "<<allCandiRegionsVec<<endl<<candiIntEng.t();
//                 cout<<candiExtEng.t();
//                 cout<<candiTotalEng.t();

                outStrm<<"Candis for "<<roiLabels[idxNets][idxROI]<<" : "<<allCandiRegionsVec<<endl<<candiIntEng.t();
                outStrm<<candiExtEng.t();
                outStrm<<candiTotalEng.t();

                allExtEng ( idxTmplRegions,idxNets ) = candiExtEng[minPox];
                allIntEng ( idxTmplRegions,idxNets ) = candiIntEng[minPox];
                allTotalEng ( idxTmplRegions,idxNets ) = candiTotalEng[minPox];
                allCorrelation ( idxTmplRegions,idxNets ) = ( 1- allExtEng ( idxTmplRegions,idxNets ) ) *sgn ( centerPrfs ( idxTmplRegions,idxNets ) );
                allIndividualFuncNetsNodes ( idxTmplRegions,idxNets ) = allCandiRegionsVec[minPox] + step;
                mapSub2tmpIdx ( allCandiRegionsVec[minPox] + step, idxNets ) = idxTmplRegions;

//                 cout<<"Results for "<<roiLabels[idxNets][idxROI]<<" : "<<right << setw ( 6 )  << idxTmplRegions << setw ( 6 ) << allCandiRegionsVec[minPox] << setw ( 10 ) << right << candiIntEng[minPox] << setw ( 10 ) << candiExtEng[minPox] << setw ( 10 ) << candiTotalEng[minPox] << endl;
                outStrm << right << setw ( 6 )  << idxTmplRegions << setw ( 6 ) << allCandiRegionsVec[minPox] << setw ( 10 ) << right << candiIntEng[minPox] << setw ( 10 ) << candiExtEng[minPox] << setw ( 10 ) << candiTotalEng[minPox] << endl;
//                 cout<<"--------------------------------------------------------------------------------------------------------------------"<<endl;
                outStrm<<"--------------------------------------------------------------------------------------------------------------------"<<endl;
            }

        }// end for loop::idxROI;
    } // end for loop::idxNets


    //output results;
    //cout<<"save results...\n";
    string outbase ( subjectRootDir+argv[4] );
    allIntEng.save ( outbase+".intEng" );
    allExtEng.save ( outbase+".extEng" );
    allTotalEng.save ( outbase+".totEng" );
    allIndividualFuncNetsNodes.save ( outbase+".indNodes" );
    allCorrelation.save ( outbase+".indProfiles" );

    //rendering surfaces for individual nets;
    if ( writeIndivFuncNetOpt.set() ) {
        CTriSurface lh ( subjectRootDir+parcelConfigFileName[1] ), rh ( subjectRootDir+parcelConfigFileName[3] );
        string attriBaseName ( "net." );
        for ( int idxNet = 0; idxNet < centerPrfs.n_cols ; ++idxNet ) {
            string attriName = attriBaseName;
            attriName+= idxNet;

            vector<float> attributeLh ( lh.GetNumOfPoints(),0 ), attributeRh ( rh.GetNumOfPoints(),0 );
            //cout<<"#################################### "<< idxNet<<" #############################"<<endl;
            outStrm<<"#################################### "<< idxNet<<" #############################"<<endl;
            //for lh;
            for ( int idxPartition = 0; idxPartition < centerPrfs.n_rows/2 ; ++idxPartition ) {
                CPartition& part = subjParts[idxPartition];
                vector<int>& pids = part.parcelIds;
                int subjROIidx = part.partitionLabel;
                int tmplROIidx = mapSub2tmpIdx ( subjROIidx,idxNet );

                if ( -1== tmplROIidx )
                    continue;
                if ( abs ( centerPrfs ( tmplROIidx,idxNet ) ) > minThreshProfileOpt.value() ) {
                    for ( int idxPid = 0; idxPid < pids.size() ; ++idxPid ) {
                        attributeLh[pids[idxPid]] = allCorrelation ( tmplROIidx,idxNet );
                    } // end for loop::idxPid
                    //cout<< subjROIidx<<" ";
                    outStrm<< subjROIidx<<" ";
                }

            } // end for loop::idxPartition

            //for rh;
            for ( int idxPartition = centerPrfs.n_rows/2; idxPartition < centerPrfs.n_rows ; ++idxPartition ) {
                CPartition& part = subjParts[idxPartition];
                vector<int>& pids = part.parcelIds;
                int subjROIidx = idxPartition;
                int tmplROIidx = mapSub2tmpIdx ( subjROIidx,idxNet );
                if ( -1== tmplROIidx )
                    continue;
                if ( abs ( centerPrfs ( tmplROIidx,idxNet ) ) > minThreshProfileOpt.value() ) {
                    for ( int idxPid = 0; idxPid < pids.size() ; ++idxPid ) {
                        attributeRh[pids[idxPid]] = allCorrelation ( tmplROIidx,idxNet );
                    } // end for loop::idxPid

                    //cout<< subjROIidx<<" ";
                    outStrm<< subjROIidx<<" ";
                }

            } // end for loop::idxPartition

            lh.AddPointDataScalar ( attributeLh,attriName );
            rh.AddPointDataScalar ( attributeRh,attriName );

            //cout<<endl;
            outStrm<<endl;
        } // end for loop::idxNet

        lh.CombineSurf ( rh );
        lh.SaveAs ( outbase+".FuncNets.vtk" );
    }

    //cout<<"done!"<<endl;
    outStrm<<"done!"<<endl;
    outStrm.close();
    return 0;
} //end of main function;

