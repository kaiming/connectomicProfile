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
#include "triSurface.h"

using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace KML;
using namespace std;
using  namespace Utilities;

string appName(" apply an affine transformation to a surface");
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


	CTriSurface mySurf(argv[1]); 
	Matrix allVertices(4,mySurf.GetNumOfPoints()); 
	allVertices=1;
	
	VectorType offsetIn, offsetOut; 
	KML::GetOffsetFromMHD(argv[3],offsetIn);
	KML::GetOffsetFromMHD(argv[4],offsetOut);
	

	for (int idxVertices = 0 ; idxVertices < mySurf.GetNumOfPoints(); ++idxVertices) {
	  allVertices(1,idxVertices+1) = mySurf.GetPointComponent(idxVertices,0);
	  allVertices(2,idxVertices+1) = mySurf.GetPointComponent(idxVertices,1);
	  allVertices(3,idxVertices+1) = mySurf.GetPointComponent(idxVertices,2);
	  allVertices(1,idxVertices+1)-= offsetIn.x;
	  allVertices(2,idxVertices+1)-= offsetIn.y;
	  allVertices(3,idxVertices+1)-= offsetIn.z;
	} // end of for::idxVertices
	
        Matrix mat(4,4);
	mat=read_ascii_matrix(argv[2]);
	
	Matrix newVertices = mat*allVertices; 
	
	for (int idxVertices = 0 ; idxVertices < mySurf.GetNumOfPoints(); ++idxVertices) {
	  
	  VectorType newCoords(newVertices(1, idxVertices+1), newVertices(2, idxVertices+1),newVertices(3, idxVertices+1));
	  newCoords+= offsetOut; 
	   mySurf.ChangePointCoords(idxVertices,newCoords);
	} // end of for::idxVertices

	mySurf.SaveAs(argv[5]);

	return 0;
} //end of main function;

