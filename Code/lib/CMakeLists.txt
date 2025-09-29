add_library(GrapheneLib STATIC "src/Constants.cpp" "src/funkcje.cpp" "src/Utilities.cpp")

target_include_directories(GrapheneLib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)