add_executable(values_Magnetism main.cpp)

target_link_libraries(values_Magnetism PRIVATE GrapheneLib)

set_target_properties(values_Magnetism PROPERTIES
    OUTPUT_NAME main
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})