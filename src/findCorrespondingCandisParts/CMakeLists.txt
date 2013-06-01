set(projectName findCorrespondingCandis)
cmake_minimum_required(VERSION 2.6)
set(AllLibs VTK MYFSL Alglib GSL Arma KML Metis)

foreach ( pkg ${AllLibs} )
find_package(${pkg} REQUIRED)
if(${pkg}_FOUND)
  include( ${USE_${pkg}_FILE} )
endif(${pkg}_FOUND)
endforeach(pkg)



PROJECT(${projectName})

IF(libName)
ADD_LIBRARY(${libName} ${LIB_TYPE} ${libName}.cpp)
ENDIF(libName)

ADD_EXECUTABLE(${projectName} ${projectName}.cpp)

TARGET_LINK_LIBRARIES(${projectName} ${libName} ${KML_LIBRARIES} ${Graclus_LIBRARIES} ${Metis_LIBRARIES} ${Arma_LIBRARIES} ${GSL_LIBRARIES} ${MYFSL_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${Alglib_LIBRARIES})
 
INSTALL(TARGETS ${projectName} DESTINATION bin)


IF(libName)
INSTALL(TARGETS ${libName} DESTINATION  lib)
FILE(GLOB allHeader "*.h") 
INSTALL(FILES ${allHeader} DESTINATION include)
ENDIF(libName)


