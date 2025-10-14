cmake_minimum_required(VERSION 3.10)
project(GrapheneCalculations LANGUAGES CXX Fortran)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

# compilation flags
if(NOT CMAKE_BUILD_TYPE)
message(FATAL_ERROR "You need to specify build type: -DCMAKE_BUILD_TYPE=<Debug/Release>")
endif()
message(STATUS "If CmakeCache already exists this check fails so be sure to provide -DCMAKE_BUILD_TYPE= everytime!")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -Wall -Wextra -DNDEBUG")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC -fPIE -warn all")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# library
add_subdirectory(lib)

# subprojects
add_subdirectory(examples/Bilayer_Magnetism_densitiesCalculation)
add_subdirectory(examples/Bilayer_noMagnetism)
add_subdirectory(examples/Bilayer_transport_v0)
