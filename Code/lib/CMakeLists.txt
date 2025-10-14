add_library(GrapheneLib STATIC 
            "src/Constants.cpp"
            "src/Bilayer.cpp"
            "src/Utilities.cpp"
            "src/Bilayer_c_interface.cpp"
            "src/Bilayer_f90_interface.f90")

set_target_properties(GrapheneLib PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mod
)

target_include_directories(GrapheneLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_include_directories(GrapheneLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mod)

