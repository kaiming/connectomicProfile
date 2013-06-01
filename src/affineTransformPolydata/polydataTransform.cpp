#include "fibers.h"
#include "iostream"
#include "options.h"
#include "newimageall.h"
#include "Vector3D.h"
#include <string>
#include <queue>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include "vtkTransformPolyDataFilter.h"
#include "vtkMatrix4x4.h"
#include "newmat.h"
#include "newmatio.h"
#include "kaimingCommon.h"

using namespace NEWMAT;
using namespace NEWIMAGE;
using namespace KML;
using namespace std;
using  namespace Utilities;

string appName(" apply an affine transformation to a polydata fibers  or surfaces");
string appExample("applicationName inputPolydata.vtk matrix.mat mhdFileIn  mhdFileRef output.vtk [optoins] ");
Option<bool> helpOpt(string("-h,--help"),false,string("display this help information, default: false( don't display"),false,requires_argument);


int main(int argc, char **argv) {

	OptionParser cmdParser(appName,appExample);
	cmdParser.add(helpOpt);
	cmdParser.parse_command_line(argc, argv, 3);
	if (6!= argc)
	{
		cmdParser.usage();
		exit(EXIT_FAILURE);
	}




	Matrix mat(4,4);
	mat=read_ascii_matrix(argv[2]);
//	cout<<mat<<std::endl;
//	Real inputValue[] = {0.999861,-0.0400237,-0.0737226,7.60398,0.0405115,0.935613,0.014767,-10.1322,0.0995461,-0.0250695,0.961496,-2.46203,0,0,0,1};
//	mat<<inputValue;

	vtkMatrix4x4* vtkMat=vtkMatrix4x4::New();
	for (int line = 0; line < 4; ++line) {
		for (int col = 0; col < 4; ++col) {
			vtkMat->SetElement(line,col, mat(line+1,col+1));
		} //end loop:: col
	} //end loop:: line

	vtkPolyDataReader*  surfReader1= vtkPolyDataReader::New();
	surfReader1->SetFileName(argv[1]);
	surfReader1->Update();

	vtkPolyData* surfData1= vtkPolyData::New();
	surfData1= surfReader1->GetOutput();
	surfData1->Update();

	Vector3D<float> offsetIn;
	GetOffsetFromMHD(argv[3], offsetIn);
	Vector3D<float> offsetRef;
	GetOffsetFromMHD(argv[4], offsetRef);
//	cout<<offsetRef<<endl;

	//remove offset before transform;
	vtkPoints* thePoints= vtkPoints::New();
	thePoints= surfData1->GetPoints();
	for (int index = 0; index < surfData1->GetNumberOfPoints(); ++index) {
		double coord[3];
		thePoints->GetPoint(index, coord);
		coord[0]-=offsetIn.x;
		coord[1]-=offsetIn.y;
		coord[2]-=offsetIn.z;
		thePoints->SetPoint(index, coord);
	} //end loop:: index

	surfData1->SetPoints(thePoints);
//	{
//		vtkPolyDataWriter* writer= vtkPolyDataWriter::New();
//		writer->SetInput(surfData1);
////		writer->SetFileTypeToBinary();
////		writer->SetFileName("noOffset.vtk");
//		writer->Update();
//	}

	vtkTransform* myTransform= vtkTransform::New();
	myTransform->SetMatrix( vtkMat);

	vtkTransformPolyDataFilter* plyTransform= vtkTransformPolyDataFilter::New();
	plyTransform->SetInput( surfData1);
	plyTransform->SetTransform(myTransform);


	vtkPolyData* outData= plyTransform->GetOutput();
	outData->Update();
	//add offset after transorm;
	{
		vtkPoints* thePoints= vtkPoints::New();
		thePoints= outData->GetPoints();
		for (int index = 0; index < surfData1->GetNumberOfPoints(); ++index) {
			double coord[3];
			thePoints->GetPoint(index, coord);
			coord[0]+=offsetRef.x;
			coord[1]+=offsetRef.y;
			coord[2]+=offsetRef.z;
			thePoints->SetPoint(index, coord);
		} //end loop:: index
		outData->SetPoints(thePoints);
	}

	vtkPolyDataWriter* writer= vtkPolyDataWriter::New();
	writer->SetInput(outData);
	writer->SetFileTypeToBinary();
	writer->SetFileName(argv[5]);
	writer->Update();


	return 0;
} //end of main function;

