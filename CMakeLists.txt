cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(Distance)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable (Distance Distance.cpp)
 
if(VTK_LIBRARIES)
  target_link_libraries(Distance ${VTK_LIBRARIES})
else()
  target_link_libraries(Distance vtkHybrid vtkWidgets)
endif()