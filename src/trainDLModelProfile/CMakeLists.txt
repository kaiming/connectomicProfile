set(projectName trainfMRIDicccol)
set(libName fmriDicccol)
cmake_minimum_required(VERSION 2.6)


set(AllLibs VTK MYFSL Alglib GSL Arma KML Spams Metis)
#set(CMAKE_CXX_FLAGS "-fpermissive")
 
foreach ( pkg ${AllLibs} )
find_package(${pkg} REQUIRED)
if(${pkg}_FOUND)
  include( ${USE_${pkg}_FILE} )
endif(${pkg}_FOUND)
endforeach(pkg)

#set(USE_BLAS_LIB)

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

PROJECT(${projectName})

IF(libName)
ADD_LIBRARY(${libName} ${LIB_TYPE} ${libName}.cpp)
ENDIF(libName)

ADD_EXECUTABLE(${projectName} ${projectName}.cpp)
TARGET_LINK_LIBRARIES(${projectName} ${libName} ${KML_LIBRARIES} ${Graclus_LIBRARIES} ${Metis_LIBRARIES} 
		      ${Arma_LIBRARIES} ${GSL_LIBRARIES} ${MYFSL_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES}
		      ${Alglib_LIBRARIES} ${Spams_LIBRARIES})
INSTALL(TARGETS ${projectName} DESTINATION bin)

IF(libName)
INSTALL(TARGETS ${libName} DESTINATION   ${KML_DIR}/lib)
FILE(GLOB allHeader "*.h") 
INSTALL(FILES ${allHeader} DESTINATION  ${KML_DIR}/include)
ENDIF(libName)



