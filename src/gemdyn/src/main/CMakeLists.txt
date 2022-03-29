message(STATUS "(EC) Generating gemdyn binaries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

if (DEFINED ENV{EC_ARCH})
  if (NOT DEFINED ENV{CRAYPE_VERSION})
    set(EXTRA_LIBRARIES "-lopen-pal")
  endif()
endif()

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gemdyn
   rpnphy
   rpnphy_cpl_stubs
   modelutils
   modelutils_tmg_stubs
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(checkdmpart gemdyn rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gemdyn
   rpnphy
   rpnphy_cpl_stubs
   modelutils
   modelutils_tmg_stubs
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(gemgrid gemdyn rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(prgemnml mainprgemnml.F90)
target_link_libraries(prgemnml
   gemdyn
   rpnphy
   rpnphy_cpl_stubs
   modelutils
   modelutils_tmg_stubs
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(prgemnml gemdyn rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)

install(TARGETS checkdmpart gemgrid prgemnml RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

