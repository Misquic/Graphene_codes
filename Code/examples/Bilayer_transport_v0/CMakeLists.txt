add_executable(Bilayer_transport_v0 ${CMAKE_CURRENT_SOURCE_DIR}/main.f90)
target_sources(Bilayer_transport_v0 PRIVATE ../TransportUtils.f90)

# all needed libraries (mainly for bubel)
target_link_libraries(Bilayer_transport_v0 PRIVATE GrapheneLib
PUBLIC bubel
#  superlu_fortran
#  superlu
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
### target_link_options(Bilayer_transport_v0 PUBLIC -lmkl_rt)
target_include_directories(Bilayer_transport_v0 PUBLIC ~/libs/include/ ~/libs/include/suitesparse ~/libs/include/Bubel/ )
#~/libs/include/FORTRAN/)
target_link_directories(Bilayer_transport_v0 PUBLIC ${MKLROOT}/lib/intel64 ~/libs/lib64)


set_target_properties(Bilayer_transport_v0 PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(Bilayer_transport_v0 PROPERTIES
  OUTPUT_NAME Transport2D
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
