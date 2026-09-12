add_executable(3D_Transport_multilead ${CMAKE_CURRENT_SOURCE_DIR}/main.f90)
target_sources(3D_Transport_multilead PRIVATE ../TransportUtils.f90)

# all needed libraries (mainly for bubel)
target_link_libraries(3D_Transport_multilead PRIVATE GrapheneLib
  PUBLIC bubel
  umfpack
  cholmod
  amd
  camd
  colamd
  ccolamd
  suitesparseconfig
  mkl_blas95_lp64
  mkl_lapack95_lp64
  mkl_rt)

# all needed libraries (for bubel)
target_include_directories(3D_Transport_multilead PUBLIC ~/libs/include/ ~/libs/include/suitesparse ~/libs/include/Bubel/ )
target_link_directories(3D_Transport_multilead PUBLIC ${MKLROOT}/lib/intel64 ~/libs/lib64)


set_target_properties(3D_Transport_multilead PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(3D_Transport_multilead PROPERTIES
  OUTPUT_NAME Transport3D
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
