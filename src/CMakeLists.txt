
cmake_minimum_required(VERSION 2.8)
set(LIB_TYPE SHARED)



find_package(Gifti REQUIRED)
if(Gifti_FOUND)
  include( ${USE_Gifti_FILE} )
endif(Gifti_FOUND)






set(subdirs affineTransformPolydata  fineParcellation   trainDLModelProfile  calculateConnectivityPattern  fsConverter  visualizeFuncNets findCommonRsns generateCentralSurface  volumeConverter findCorrespondingCandisParts  getIndividualProfile findFsRegMatchingVertices  rmCentralSurfSpikes fineParcellationHCP)

foreach ( program ${subdirs} )
add_subdirectory (${program} ) 
endforeach(program)



