add_executable(Bilayer_noMagnetism main.cpp)

target_link_libraries(Bilayer_noMagnetism PRIVATE GrapheneLib)

set_target_properties(Bilayer_noMagnetism PROPERTIES
    OUTPUT_NAME main
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})