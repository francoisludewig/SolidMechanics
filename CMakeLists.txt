cmake_minimum_required(VERSION 3.6)

project(SolidMechanicsProgram VERSION 0.1.0)

add_definitions(-std=c++11) 


#LINK_DIRECTORIES(${LIB_DIR}/usr/local/lib/../../)

set(SOURCES_FILES
	main.cpp
)

set(FILES
    ${SOURCES_FILES}
)

add_executable(SolidMechanicsProgram
	${FILES}
)

target_link_libraries(
	SolidMechanicsProgram
	GeometricalSpaceObjects.libs	
)

enable_testing()

add_subdirectory("GeometricalSpaceObjects")
add_subdirectory("GeometricalSolid")
