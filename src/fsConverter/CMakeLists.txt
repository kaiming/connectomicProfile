set(projectName "fsutils")

set(AllLibs VTK ITK MYFSL Alglib GSL Arma KML Metis)

foreach ( pkg ${AllLibs} )
find_package(${pkg} REQUIRED)
if(${pkg}_FOUND)
  include( ${USE_${pkg}_FILE} )
endif(${pkg}_FOUND)
endforeach(pkg)


cmake_minimum_required(VERSION 2.6) #cmake version requirement; 


#set(CMAKE_CXX_FLAGS "-fopenmp -lgomp")
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
 
LINK_LIBRARIES(${KML_LIBRARIES} ${Graclus_LIBRARIES} ${Metis_LIBRARIES} ${Arma_LIBRARIES} ${GSL_LIBRARIES} 
		      ${MYFSL_LIBRARIES} ${VTK_LIBRARIES} ${ITK_LIBRARIES} ${Alglib_LIBRARIES} ${MATLAB_LIBRARIES} )
		      
ADD_EXECUTABLE(fsRas2VoxAddOri fsRas2VoxAddOri.cpp)
 
ADD_EXECUTABLE(fsRas2Vox fsRas2Vox.cpp)
 
ADD_EXECUTABLE(fs2vtk fs2vtk.cpp)
 
ADD_EXECUTABLE(fsVox2Ras fsVox2Ras.cpp)
 
ADD_EXECUTABLE(fsVox2RasRmOri fsVox2RasRmOri.cpp)
 
ADD_EXECUTABLE(fsCombineLeftRightSurfs fsCombineLeftRightSurfs.cpp)
 
add_executable(fsGeneraterCentralSurf generateCS.cpp)
 
ADD_EXECUTABLE(fsCombineLeftRightCentralSurfs fsCombineLeftRightCentralSurfs.cpp)
  

 
install(TARGETS fsGeneraterCentralSurf fs2vtk fsRas2VoxAddOri fsRas2Vox fsVox2Ras fsVox2RasRmOri fsCombineLeftRightSurfs fsCombineLeftRightCentralSurfs DESTINATION  bin)
 
 


