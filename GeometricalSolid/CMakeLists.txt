cmake_minimum_required(VERSION 3.6)

set(HEADER_FILES
  Include/Solid.h
  Include/Shape.h
  Include/Sphere.h
  Include/Rectangle.h
  Include/Disk.h
  Include/Cone.h
  Include/Elbow.h
  Include/Formatter/SolidFormatter.h
  Include/Parser/SolidParser.h
)

set(SOURCE_FILES
  Source/Solid.cpp
)

add_library(GeometricalSolid.libs
 ${SOURCE_FILES}
 ${HEADER_FILES}
)

target_link_libraries(
  GeometricalSolid.libs	
  GeometricalSpaceObjects.libs
)


target_include_directories(
GeometricalSolid.libs INTERFACE
Include
Include/Formatter
Include/Parser
)



add_subdirectory("Tests")
