add_executable(Bilayer_Magnetism_densities main.cpp)

target_link_libraries(Bilayer_Magnetism_densities PRIVATE GrapheneLib)

set_target_properties(Bilayer_Magnetism_densities PROPERTIES
    OUTPUT_NAME main
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

