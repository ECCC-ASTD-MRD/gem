message(STATUS "(EC) Generating maingem* Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILER_FLAGS}")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILER_FLAGS}")
endif()

add_executable(maingemdm maingemdm.F90)

target_compile_definitions(maingemdm PUBLIC
  dyn_VERSION="${dyn_VERSION}"
  modelutils_VERSION="${modelutils_VERSION}"
  phy_VERSION="${phy_VERSION}"
)

target_link_libraries(maingemdm
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   modelutils_ov_ifort
   vgrid::vgrid
   tdpack::tdpack
   rpncomm::rpncomm
   rmn::rmn-ompi
   OpenMP::OpenMP_Fortran
   MPI::MPI_Fortran
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

install(TARGETS maingemdm RUNTIME DESTINATION bin)
