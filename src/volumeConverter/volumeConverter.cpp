#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExceptionObject.h"

# include "utils/options.h"
# include <fstream>
# include <iostream>
# include <string>
# include <cstdlib>
# include <cstdio>
# include <cassert>


using namespace std;
using namespace Utilities;
using namespace itk;

 


int main(int argc, char** argv)
{

 
  
	string appName("Convert volume between different formats");
	string example("volumeConverter <inputFile.ext> <outputFile.ext> [-t voxelType]");

	Option<bool> help(string("-h,--help"),false,string("displays this help, then exits"),false,no_argument);
	Option<string> typeOpt(string("-t,--type"),"short",string("the data type, valid types: char, uchar, short, int,float,double, default: short"),false,requires_argument);
	Option<int> dimOpt("-d",3,string("the dimension of the data, default: 3"), false,requires_argument);

	OptionParser cmdOptions(appName,example);
	cmdOptions.add(typeOpt);
	cmdOptions.add(help);
 	cmdOptions.add(dimOpt);
 


	if (argc<3)
	{
		if (1==argc)
		{
			cmdOptions.usage();
			exit(EXIT_FAILURE);
		}
		string cmdLine(argv[1]); 
		if (cmdLine.find("-h")!=string::npos || cmdLine.find("--help")!=string::npos)
		{
			cmdOptions.usage();
			exit(EXIT_FAILURE);
		}
		cerr<<"error: not enough arguments, use -h for help "<<endl;
		exit(EXIT_FAILURE);
	}

	cmdOptions.parse_command_line(argc,argv,2);

	//const unsigned int dim( dimOpt.value());
	const unsigned int dim=3; 
	if (typeOpt.value()== string("char"))
	{
		typedef char VoxelType;
		typedef itk::Image<VoxelType,dim> InputVolType;
		typedef itk::Image<VoxelType,dim> OutputVolType;

		typedef itk::ImageFileReader<InputVolType> VolReaderType;
		typedef itk::ImageFileWriter<OutputVolType> VolWriterType; 

		VolReaderType::Pointer reader= VolReaderType::New();
		VolWriterType::Pointer  writer= VolWriterType::New();


		reader->SetFileName(argv[1]);
		writer->SetFileName(argv[2]);
		writer->SetInput(reader->GetOutput());

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}
	else if ( typeOpt.value()==string("uchar"))
	{
		typedef unsigned char VoxelType;
		typedef itk::Image<VoxelType,dim> InputVolType;
		typedef itk::Image<VoxelType,dim> OutputVolType;

		typedef itk::ImageFileReader<InputVolType> VolReaderType;
		typedef itk::ImageFileWriter<OutputVolType> VolWriterType; 

		VolReaderType::Pointer reader= VolReaderType::New();
		VolWriterType::Pointer  writer= VolWriterType::New();


		reader->SetFileName(argv[1]);
		writer->SetFileName(argv[2]);
		writer->SetInput(reader->GetOutput());

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}
	else if (typeOpt.value()==string("short"))
	{
		typedef short VoxelType;
		typedef itk::Image<VoxelType,dim> InputVolType;
		typedef itk::Image<VoxelType,dim> OutputVolType;

		typedef itk::ImageFileReader<InputVolType> VolReaderType;
		typedef itk::ImageFileWriter<OutputVolType> VolWriterType; 

		VolReaderType::Pointer reader= VolReaderType::New();
		VolWriterType::Pointer  writer= VolWriterType::New();


		reader->SetFileName(argv[1]);
		writer->SetFileName(argv[2]);
		writer->SetInput(reader->GetOutput());

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}
	else if (typeOpt.value()==string("int"))
	{
		typedef int VoxelType;
		typedef itk::Image<VoxelType,dim> InputVolType;
		typedef itk::Image<VoxelType,dim> OutputVolType;

		typedef itk::ImageFileReader<InputVolType> VolReaderType;
		typedef itk::ImageFileWriter<OutputVolType> VolWriterType; 

		VolReaderType::Pointer reader= VolReaderType::New();
		VolWriterType::Pointer  writer= VolWriterType::New();


		reader->SetFileName(argv[1]);
		writer->SetFileName(argv[2]);
		writer->SetInput(reader->GetOutput());

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}
	else if (typeOpt.value()==string("float"))
	{
		typedef float VoxelType;
		typedef itk::Image<VoxelType,dim> InputVolType;
		typedef itk::Image<VoxelType,dim> OutputVolType;

		typedef itk::ImageFileReader<InputVolType> VolReaderType;
		typedef itk::ImageFileWriter<OutputVolType> VolWriterType; 

		VolReaderType::Pointer reader= VolReaderType::New();
		VolWriterType::Pointer  writer= VolWriterType::New();


		reader->SetFileName(argv[1]);
		writer->SetFileName(argv[2]);
		writer->SetInput(reader->GetOutput());

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}
	else if (typeOpt.value()==string("double"))
	{
		typedef double VoxelType;
		typedef itk::Image<VoxelType,dim> InputVolType;
		typedef itk::Image<VoxelType,dim> OutputVolType;

		typedef itk::ImageFileReader<InputVolType> VolReaderType;
		typedef itk::ImageFileWriter<OutputVolType> VolWriterType; 

		VolReaderType::Pointer reader= VolReaderType::New();
		VolWriterType::Pointer  writer= VolWriterType::New();


		reader->SetFileName(argv[1]);
		writer->SetFileName(argv[2]);
		writer->SetInput(reader->GetOutput());

		try
		{
			writer->Update();
		}
		catch( itk::ExceptionObject & err )
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}
	else
	{
		cout<<"Invalid data type!"<<endl;
		cmdOptions.usage();
		exit(EXIT_FAILURE);
	}



	
	return 0;
}
