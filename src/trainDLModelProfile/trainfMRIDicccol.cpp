/*
*	copyright by kaiming li (kaiming.li@emory.edu) 2012.
*	License : GPL v3
*
*/
#include <string>
#include <iostream>
#include "options.h"
#include "kaimingCommon.h"
#include "fmriDicccol.h"
#include "apClustering.h"

using namespace arma;


using namespace Utilities;
using namespace KML;
using namespace std;


int main ( int argc, char **argv )
{

        string appName ( argv[0] );
        string appExample ( "partitionConfigFile outputBaseName [ options] \
			    \n\nDetails:\n\
			    \tthis program reads a series of connectivity pattern files, and trains a connectomic profile for each ROI\n\
			    \tpartitionConfigFile contains the names of output from the pipeline.\n\
			    \tgrpProfile is the output for resultant profiles, written in arma format\n" );
        Option<bool> helpOpt ( string ( "-h" ),false,string ( "display this help information. " ),false,no_argument );
        Option<int> modeOpt ( "--mode",3,"the mode of engergy functions, e.g., L1COEFFS(0), L2ERROR(1), PENALTY(2), SPARSITY(3), L2ERROR2(4), PENALTY2(5), default: 3", false, requires_argument );
        Option<int> modeDOpt ( "--modeD",3,"the constraints on dictionary:  L2(0),  L1L2(1), L1L2FL(2), L1L2MU(3), default: 3",false,requires_argument );
        Option<int> modePOpt ( "--modeP",0,"the method to minimization:  AUTO(0), PARAM1(1), PARAM2(2), PARAM3(3), default: 0",false,requires_argument );
        Option<float> lambdaOpt ( "--lamda",1," lambda, default: 1", false,requires_argument );
        Option<float> lambda2Opt ( "--lamda2", 0.15, "lambda2, default: 0", false,requires_argument );
        Option<float> gamma1Opt ( "--gamma1",0.4,"gamma1, default:0.4",false,requires_argument );
        Option<float> gamma2Opt ( "--gamma2",0,"gamma2, default:0",false,requires_argument );
        Option<bool> posAlphaOpt ( "--posAlpha",true,"where or not constrain alpha to be positive, default: true",false,no_argument );
        Option<bool> posDOpt ( "--posD",false,"whether or not constrain Dictionary to be positive, default: false",false,no_argument );
        Option<int> iterOpt ( "--iter", 300,"the number of iteration, default: 300",false,requires_argument );
        Option<bool> verboseOpt ( "-v",false,"verbose or not,default, false",false,no_argument );
        Option<int> numDicOpt ( "-K",5,"the number of dic elems, default: 5",false,requires_argument );
        Option<int> batchSizeOpt ( "-b",30,"the number of batch, default: 30", false,requires_argument );
        Option<string> aplibFullNameOpt ( "--aplib","/home/kli/bin/apclusterunix64.so","the full path name of ap lib;", false,requires_argument );
        Option<float> prefrenceOpt ( "-p",0,"the prefrence of ap clustering, if not set, will use mean similarity.", false,requires_argument );
        Option<string> leftPartCenterSphereOpt ( "--leftParts","lh.resx.200.CentersonSphere.vtk","the name of left part centers at reg space, default: lh.parts.200.CentersonSphere.vtk",false,requires_argument );
        Option<string> rightPartCenterSphereOpt ( "--rightParts","rh.resx.200.CentersonSphere.vtk","the name of right part centers at reg space, default: rh.parts.200.CentersonSphere.vtk",false,requires_argument );
        Option<string> leftLabelMatchAtSubjOpt ( "--leftMatch","labelMatch.lh.resx.200.tmplAtSubj", "the name of left label match for tmpl at subject space, default: labelMatch.lh.resx.200.tmplAtSubj", false,requires_argument );
        Option<string> rightLabelMatchAtSubjOpt ( "--rightMatch","labelMatch.rh.resx.200.tmplAtSubj", "the name of right label match for tmpl at subject space, default: labelMatch.rh.resx.200.tmplAtSubj", false,requires_argument );
        Option<string> subjDirsListOpt ( "-L","","the list of allReg directory for subjects, no default, must be set", true, requires_argument );

        OptionParser cmdParser ( appName,appName+"\t"+appExample );
        cmdParser.add ( helpOpt );
        cmdParser.add ( verboseOpt );
        cmdParser.add ( numDicOpt );
        cmdParser.add ( batchSizeOpt );
        cmdParser.add ( prefrenceOpt );
        cmdParser.add ( modeOpt );
        cmdParser.add ( modePOpt );
        cmdParser.add ( modeDOpt );
        cmdParser.add ( lambdaOpt );
        cmdParser.add ( lambda2Opt );
        cmdParser.add ( gamma2Opt );
        cmdParser.add ( gamma1Opt );
        cmdParser.add ( posDOpt );
        cmdParser.add ( posAlphaOpt );
        cmdParser.add ( iterOpt );
        cmdParser.add ( aplibFullNameOpt );
        cmdParser.add ( leftLabelMatchAtSubjOpt );
        cmdParser.add ( rightLabelMatchAtSubjOpt );
        cmdParser.add ( leftPartCenterSphereOpt );
        cmdParser.add ( rightPartCenterSphereOpt );
        cmdParser.add ( subjDirsListOpt );

        cmdParser.parse_command_line ( argc,argv,2 );

        if ( 3 > argc || cmdParser.check_compulsory_arguments ( true ) == false ) {
                cmdParser.usage();
                exit ( EXIT_FAILURE );
        }


        //set up parameters for dic learning;
        ParamDictLearn<float> parameters;
        parameters.mode = SPAMS::constraint_type ( modeOpt.value() );
        parameters.posAlpha = posAlphaOpt.value();
        parameters.lambda= lambdaOpt.value();
        parameters.modeD=SPAMS::constraint_type_D ( modeDOpt.value() );
        parameters.gamma1=gamma1Opt.value();
        parameters.modeParam=SPAMS::mode_compute ( modePOpt.value() );
        parameters.gamma2 = gamma2Opt.value(); // for modeD=2;
        parameters.iter=iterOpt.value();
        parameters.verbose=verboseOpt.value();
	
	
	//setup input connectivity files; 
	vector<string> allConfigNames; 
	KML::ReadNameList(argv[1],allConfigNames);
 

        CfMRIDicccol objfMRIDicccol;
	objfMRIDicccol.SetConnectivityFileName(allConfigNames[8]);
	  objfMRIDicccol.SetNames4TmplLabelMatchAtSubj ( allConfigNames[6],allConfigNames[7]);
        objfMRIDicccol.SetNames4PartitionCenterSurfAtRegSpace ( allConfigNames[10],allConfigNames[11] );
	
        objfMRIDicccol.SetData ( subjDirsListOpt.value() );
        objfMRIDicccol.SetParamDicLearn ( parameters );
        objfMRIDicccol.SetTrainer ( numDicOpt.value(),batchSizeOpt.value(),-1 );
		
        objfMRIDicccol.LearnProfileAndCoordModels();
        objfMRIDicccol.SaveProfileAndCoordModels ( argv[2] );


        //now clustering profiles;
        list<int> v1,v2;
        list<float> sims;
        float meanSim = objfMRIDicccol.ComputeSims4Profiles ( v1,v2,sims,string ( argv[2] ) +".cluster.sims" );
        CAPClustering apc;
        apc.SetNPoints ( objfMRIDicccol.GetNROIs() );
        apc.SetAplibName ( aplibFullNameOpt.value() );
        apc.SetMaxIteration ( 4000 );
        apc.SetPreference ( prefrenceOpt.set() ? prefrenceOpt.value() : meanSim );
        apc.Clustering ( v1,v2,sims );
        std::vector< int > idx=apc.GetClusterResultIdx();
        std::vector< int > centers=apc.GetCenters();
        ivec armaIdx ( &idx[0],idx.size() );
        armaIdx.save ( string ( argv[2] ) +".cluster.idx" );


        return 0;
} //end of main function;


