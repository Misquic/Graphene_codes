add_executable(Bilayer_transport_v0 ${CMAKE_CURRENT_SOURCE_DIR}/main.f90)

# all needed libraries (mainly for bubel)
target_link_libraries(Bilayer_transport_v0 PRIVATE GrapheneLib 
  PUBLIC bubel
  mkl_lapack95_lp64
  mkl_blas95_lp64
  mkl_intel_lp64
  superlu_fortran
  superlu)
# all needed libraries (for bubel)
target_link_options(Bilayer_transport_v0 PUBLIC -qmkl)
target_include_directories(Bilayer_transport_v0 PUBLIC /usr/local/include/bubel/)
target_link_directories(Bilayer_transport_v0 PUBLIC /opt/intel/mkl/oneapi/mkl/latest/lib)


set_target_properties(Bilayer_transport_v0 PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(Bilayer_transport_v0 PROPERTIES
  OUTPUT_NAME main
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
