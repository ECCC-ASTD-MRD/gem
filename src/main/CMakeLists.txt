message(STATUS "(EC) Generating gemdyn binaries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILER_FLAGS}")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILER_FLAGS}")
endif()

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
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

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
   vgrid::vgrid
   tdpack::tdpack
   rpncomm::rpncomm
   rmn::rmn-ompi
   OpenMP::OpenMP_Fortran
   MPI::MPI_Fortran
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

add_executable(prgemnml mainprgemnml.F90)
target_link_libraries(prgemnml
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   modelutils_tmg_stubs
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

install(TARGETS checkdmpart gemgrid prgemnml RUNTIME DESTINATION bin)
