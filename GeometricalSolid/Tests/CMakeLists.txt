cmake_minimum_required(VERSION 3.1.2)

set(SOURCES_FILES
	main.cpp
  TestSolid.cpp
  TestSphere.cpp
)

set(FILES
    ${SOURCES_FILES}
)

add_executable(
	GeometricalSolid.Tests 
	${FILES}
	${tests_files}
)

find_package(GTest REQUIRED)

enable_testing()

include_directories(
    ${GTEST_INCLUDE_DIRS}
)

target_link_libraries(
	GeometricalSolid.Tests 
	GeometricalSolid.libs
	gmp
	mpfr
	${GTEST_LIBRARIES}
)

link_directories(/usr/local/lib)
