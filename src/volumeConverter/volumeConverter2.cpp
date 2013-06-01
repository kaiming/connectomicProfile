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
#include <string>
#include <iostream>
#include "options.h"
#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include <itkImageIOBase.h>
#include <itkNiftiImageIO.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include "itkCastImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"

using namespace Utilities;
using namespace std;
using namespace itk;

 
void GetImageType ( std::string fileName,ImageIOBase::IOPixelType& pixelType,ImageIOBase::IOComponentType& componentType )
{
        typedef itk::Image<unsigned char, 4> ImageType;
        itk::ImageFileReader<ImageType>::Pointer imageReader = itk::ImageFileReader<ImageType>::New();
        imageReader->SetFileName ( fileName.c_str() );
        imageReader->UpdateOutputInformation();

        pixelType = imageReader->GetImageIO()->GetPixelType();
        componentType = imageReader->GetImageIO()->GetComponentType();
}


template<typename PixelType>
int DoIt ( int argc, char* argv[], PixelType dummy)
{
        typedef PixelType VoxelType; 
        const unsigned int dim = 3; 
        typedef itk::Image<VoxelType,dim> InputVolType;
        typedef itk::Image<VoxelType,dim> OutputVolType;

        typedef itk::ImageFileReader<InputVolType> VolReaderType;
        typedef itk::ImageFileWriter<OutputVolType> VolWriterType;

        typename VolReaderType::Pointer reader= VolReaderType::New();
        typename VolWriterType::Pointer  writer= VolWriterType::New();


        reader->SetFileName ( argv[1] );
        writer->SetFileName ( argv[2] );
        writer->SetInput ( reader->GetOutput() );

        try {
                writer->Update();
        } catch ( itk::ExceptionObject & excep ) {
                std::cerr << "Exception catched !" << std::endl;
                std::cerr << excep << std::endl;
        }

        return EXIT_SUCCESS;

}
 

int main ( int argc, char * argv[] )
{

        string appName("Convert volume between different formats");
        string example("volumeConverter <inputFile.ext> <outputFile.ext> ");
        Option<bool> help(string("-h,--help"),false,string("displays this help, then exits"),false,no_argument);
        OptionParser cmdOptions(appName,example);
        cmdOptions.add(help);
         
        if ( argc<3 ) {
                if ( 1==argc ) {
                        cmdOptions.usage();
                        exit ( EXIT_FAILURE );
                }
                string cmdLine ( argv[1] );
                if ( cmdLine.find ( "-h" ) !=string::npos || cmdLine.find ( "--help" ) !=string::npos ) {
                        cmdOptions.usage();
                        exit ( EXIT_FAILURE );
                }
                cerr<<"error: not enough arguments, use -h for help "<<endl;
                exit ( EXIT_FAILURE );
        }

        cmdOptions.parse_command_line ( argc,argv,2 );

        itk::ImageIOBase::IOPixelType pixelType;
        itk::ImageIOBase::IOComponentType componentType;

        try {

                GetImageType ( std::string ( argv[1] ), pixelType, componentType );


 
                switch ( componentType ) {
#ifndef WIN32
                case itk::ImageIOBase::UCHAR:
                        std::cout<<"Input File Type: UCHAR"<<std::endl;
                        return DoIt ( argc, argv, static_cast<unsigned char> ( 0 ) );
                        break;
                case itk::ImageIOBase::CHAR:
                        std::cout<<"Input File Type: CHAR"<<std::endl;
                        return DoIt ( argc, argv, static_cast<char> ( 0 ) );
                        break;
#endif
                case itk::ImageIOBase::USHORT:
                        std::cout<<"Input File Type: USHORT"<<std::endl;
                        return DoIt ( argc, argv, static_cast<unsigned short> ( 0 ) );
                        break;
                case itk::ImageIOBase::SHORT:
                        std::cout<<"Input File Type: SHORT"<<std::endl;
                        return DoIt ( argc, argv, static_cast<short> ( 0 ) );
                        break;
                case itk::ImageIOBase::UINT:
                        std::cout<<"Input File Type: UINT"<<std::endl;
                        return DoIt ( argc, argv, static_cast<unsigned int> ( 0 ) );
                        break;
                case itk::ImageIOBase::INT:
                        std::cout<<"Input File Type: INT"<<std::endl;
                        return DoIt ( argc, argv, static_cast<int> ( 0 ) );
                        break;
#ifndef WIN32
                case itk::ImageIOBase::ULONG:
                        std::cout<<"Input File Type: ULONG"<<std::endl;
                        return DoIt ( argc, argv, static_cast<unsigned long> ( 0 ) );
                        break;
                case itk::ImageIOBase::LONG:
                        std::cout<<"Input File Type: LONG"<<std::endl;
                        return DoIt ( argc, argv, static_cast<long> ( 0 ) );
                        break;
#endif
                case itk::ImageIOBase::FLOAT:
                        std::cout<<"Input File Type: FLOAT"<<std::endl;
                        return DoIt ( argc, argv, static_cast<float> ( 0 ) );
                        break;
                case itk::ImageIOBase::DOUBLE:
                        std::cout<<"Input File Type: DOUBLE"<<std::endl;
                        return DoIt ( argc, argv, static_cast<double> ( 0 ) );
                        break;
                case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
                default:
                        std::cout << "unknown component type" << std::endl;
                        break;
                }

        } catch ( itk::ExceptionObject &excep ) {
                std::cerr << argv[0] << ": exception caught !" << std::endl;
                std::cerr << excep << std::endl;
                return EXIT_FAILURE;
        }
 
        return EXIT_SUCCESS;

}


 