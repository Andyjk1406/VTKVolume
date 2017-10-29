cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(VTKVolume)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable (VTKVolume VTKVolume.cpp)
 
if(VTK_LIBRARIES)
  target_link_libraries(VTKVolume ${VTK_LIBRARIES})
else()
  target_link_libraries(VTKVolume vtkHybrid vtkWidgets)
endif()
