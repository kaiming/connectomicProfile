
cmake_minimum_required(VERSION 2.6)
 
set(AllLibs VTK MYFSL Alglib GSL Arma Metis KML  Metis)

foreach ( pkg ${AllLibs} )
find_package(${pkg} REQUIRED)
if(${pkg}_FOUND)
  include( ${USE_${pkg}_FILE} )
endif(${pkg}_FOUND)
endforeach(pkg)



IF(${LIB_TYPE} STREQUAL "SHARED")
	SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
	SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
	SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
	"${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
	IF("${isSystemDir}" STREQUAL "-1")
		SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
	ENDIF("${isSystemDir}" STREQUAL "-1")
ENDIF(${LIB_TYPE} STREQUAL "SHARED")


#remove the following line for non_debug mode; 
INCLUDE_DIRECTORIES("./")
# SET(CMAKE_CXX_FLAGS "-std=c++0x")
# SET(CMAKE_CXX_FLAGS "-std=c++0x")   

set(projectName "fineParcellationMetis")

PROJECT(${projectName})

ADD_LIBRARY(kmlPartitions ${LIB_TYPE} kmlPartitions.cpp)
ADD_EXECUTABLE(${projectName} fineParcellation.cpp fineParcellator.cpp)
TARGET_LINK_LIBRARIES(${projectName} kmlPartitions ${KML_LIBRARIES} ${Graclus_LIBRARIES} ${Metis_LIBRARIES} ${Arma_LIBRARIES} ${GSL_LIBRARIES} ${MYFSL_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${Alglib_LIBRARIES})
 

install(TARGETS ${projectName} DESTINATION  bin)
IF(libName)
INSTALL(TARGETS ${libName} DESTINATION  lib)
FILE(GLOB allHeader "*.h") 
INSTALL(FILES ${allHeader} DESTINATION include)
ENDIF(libName)


