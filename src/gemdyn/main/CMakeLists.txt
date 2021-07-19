message(STATUS "Starting generation of gemdyn main Makefile")

set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

add_executable(maingemdm maingemdm.F90)
target_link_libraries(maingemdm
   gem
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   open-pal
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(maingemdm gem phy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gem
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   open-pal
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(checkdmpart gem phy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gem
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   open-pal
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(gemgrid gem phy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(prgemnml mainprgemnml.F90)
target_link_libraries(prgemnml
   gem
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   open-pal
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(prgemnml gem phy modelutils vgrid mpi-extras rpncomm rmn massv)

install(TARGETS maingemdm checkdmpart gemgrid prgemnml RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

