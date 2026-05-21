add_executable(Bilayer_Magnetism_voltages main.cpp)

target_link_libraries(Bilayer_Magnetism_voltages PRIVATE GrapheneLib)

set_target_properties(Bilayer_Magnetism_voltages PROPERTIES
    OUTPUT_NAME main
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})