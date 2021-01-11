message(STATUS "Starting generation of maingem Makefile")

set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

file(READ "${ROOT}/VERSION" GEM_VERSION)
string(REGEX REPLACE "\n$" "" GEM_VERSION "${GEM_VERSION}")
message(STATUS "GEM_VERSION=${GEM_VERSION}")

add_executable(maingemdm maingemdm.F90)
target_link_libraries(maingemdm
   gem
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(maingemdm gem phy modelutils vgrid rpncomm rmn massv)

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gem
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(checkdmpart gem phy modelutils vgrid rpncomm rmn massv)

add_executable(prgemnml mainprgemnml.F90)
target_link_libraries(prgemnml
   gem
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(prgemnml gem phy modelutils vgrid rpncomm rmn massv)

install(TARGETS prgemnml RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gem
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(gemgrid gem phy modelutils vgrid rpncomm rmn massv)

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(feseri phy modelutils vgrid rpncomm rmn massv)

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   phy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(prphynml phy modelutils vgrid rpncomm rmn massv)

install(TARGETS maingemdm gemgrid checkdmpart feseri prphynml RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)


