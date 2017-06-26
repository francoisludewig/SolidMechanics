cmake_minimum_required(VERSION 3.6)

add_library(GeometricalSpaceObjects.libs INTERFACE)

set(HEADER_FILES
  Include/Basis.h
  Include/Matrix.h
  Include/Quaternion.h
  Include/Vector.h
  Include/Point.h
  Include/VectorsQuaternionConverter.h

  Include/Formatter/BasisFormatter.h	
  Include/Formatter/QuaternionFormatter.h
  Include/Formatter/MatrixFormatter.h	
  Include/Formatter/VectorFormatter.h
  Include/Formatter/PointFormatter.h	

  Include/Parser/BasisParser.h		
  Include/Parser/QuaternionParser.h
  Include/Parser/MatrixParser.h		
  Include/Parser/VectorParser.h
  Include/Parser/PointParser.h
)

include_directories( GeometricalSpaceObjects.libs/include )



target_include_directories(
 GeometricalSpaceObjects.libs INTERFACE
 Include
)

target_link_libraries(
  GeometricalSpaceObjects.libs INTERFACE	
)

add_subdirectory("Tests")

add_custom_target(GeometricalSpaceObjectsDir SOURCES ${HEADER_FILES})
