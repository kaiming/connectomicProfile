#include "fibers.h"
#include "iostream"
#include "options.h"
#include "newimageall.h"
#include "Vector3D.h"
#include <string>
#include <queue> 
#include "newmat.h"
#include "newmatio.h"
#include "kaimingCommon.h"
#include "fibers.h"

using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace KML;
using namespace std;
using  namespace Utilities;



string appName(" apply an affine transformation to a polydata fibers ");
string appExample("applicationName inputPolydata.vtk matrix.mat mhdFileIn  mhdFileRef output.vtk [optoins] ");
Option<bool> helpOpt(string("-h,--help"),false,string("display this help , default: false( don't display"),false,requires_argument);


int main(int argc, char **argv) {

	OptionParser cmdParser(appName,appExample);
	cmdParser.add(helpOpt);
	cmdParser.parse_command_line(argc, argv, 3);
	if (6!= argc)
	{
		cmdParser.usage();
		exit(EXIT_FAILURE);
	}


	CFibers myFibers(argv[1]); 
	Matrix allVertices(4,myFibers.GetNumPoints()); 
	allVertices=1;
	
	VectorType offsetIn, offsetOut; 
	KML::GetOffsetFromMHD(argv[3],offsetIn);
	KML::GetOffsetFromMHD(argv[4],offsetOut);
	

	for (int idxVertices = 0 ; idxVertices < myFibers.GetNumPoints(); ++idxVertices) {
	  allVertices(1,idxVertices+1) = myFibers.GetPointCoords(idxVertices).x;
	  allVertices(2,idxVertices+1) = myFibers.GetPointCoords(idxVertices).y; 
	  allVertices(3,idxVertices+1) = myFibers.GetPointCoords(idxVertices).z; 
	  allVertices(1,idxVertices+1)-= offsetIn.x;
	  allVertices(2,idxVertices+1)-= offsetIn.y;
	  allVertices(3,idxVertices+1)-= offsetIn.z;
	} // end of for::idxVertices
	
        Matrix mat(4,4);
	mat=read_ascii_matrix(argv[2]);
	
	Matrix newVertices = mat*allVertices; 
	
	for (int idxVertices = 0 ; idxVertices < myFibers.GetNumPoints(); ++idxVertices) {
	  
	  VectorType newCoords(newVertices(1, idxVertices+1), newVertices(2, idxVertices+1),newVertices(3, idxVertices+1));
	  newCoords+= offsetOut; 
	   myFibers.SetPointCoords(idxVertices,newCoords);
	} // end of for::idxVertices

	myFibers.SaveAllFibers(argv[5]);

	return 0;
} //end of main function;

