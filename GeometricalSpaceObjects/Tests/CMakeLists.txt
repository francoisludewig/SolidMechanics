cmake_minimum_required(VERSION 3.1.2)

set(SOURCES_FILES
	main.cpp
	TestVector.cpp
	TestQuaternion.cpp
	TestBasis.cpp
	TestMatrix.cpp
	TestPoint.cpp	
)

set(FILES
    ${SOURCES_FILES}
)

add_executable(
	GeometricalSpaceObjects.Tests 
	${FILES}
	${tests_files}
)

find_package(GTest REQUIRED)

enable_testing()

include_directories(
    ${GTEST_INCLUDE_DIRS}
)

target_link_libraries(
	GeometricalSpaceObjects.Tests 
	GeometricalSpaceObjects.libs
  gmp
  mpfr
	${GTEST_LIBRARIES}
)

link_directories(/usr/local/lib)
