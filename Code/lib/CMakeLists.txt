cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(GrapheneLib STATIC
            "src/Constants.cpp"
            "src/Bilayer.cpp"
            "src/BilayerMagnetic.cpp"
            "src/Utilities.cpp"
            "src/Bilayer_c_interface.cpp"
            "src/Bilayer_f90_interface.f90")

set_target_properties(GrapheneLib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mod
)

target_link_libraries(GrapheneLib PUBLIC stdc++exp)

target_include_directories(GrapheneLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(GrapheneLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mod)

