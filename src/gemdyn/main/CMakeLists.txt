message(STATUS "(EC) Generating maingem Makefile")

set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

# Adding specific flags for GEM
#include(${CMAKE_SOURCE_DIR}/src/flags-GEM.cmake)

message(STATUS "(EC) CMAKE_C_FLAGS for GEM: ${CMAKE_C_FLAGS}")
message(STATUS "(EC) CMAKE_Fortran_FLAGS for GEM: ${CMAKE_Fortran_FLAGS}")

add_executable(maingemdm maingemdm.F90)
target_link_libraries(maingemdm
   gem
   phy
   modelutils
   vgrid
   rpn_comm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_executable(checkdmpart maincheckdmpart.F90)
target_link_libraries(checkdmpart
   gem
   phy
   modelutils
   vgrid
   rpn_comm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_executable(gemgrid maingemgrid.F90)
target_link_libraries(gemgrid
   gem
   phy
   modelutils
   vgrid
   rpn_comm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_executable(prgemnml mainprgemnml.F90)
target_link_libraries(prgemnml
   gem
   phy
   modelutils
   vgrid
   rpn_comm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_dependencies(prgemnml gem phy modelutils ${VGRID_LIBRARY_SRC} ${MPIEXTRAS_LIBRARY_SRC} ${RPNCOMM_LIBRARY_SRC} ${RMN_LIBRARY_SRC} ${MASSV_LIBRARY_SRC})
add_dependencies(gemgrid gem phy modelutils ${VGRID_LIBRARY_SRC} ${MPIEXTRAS_LIBRARY_SRC} ${RPNCOMM_LIBRARY_SRC} ${RMN_LIBRARY_SRC} ${MASSV_LIBRARY_SRC})
add_dependencies(maingemdm gem phy modelutils ${VGRID_LIBRARY_SRC} ${MPIEXTRAS_LIBRARY_SRC} ${RPNCOMM_LIBRARY_SRC} ${RMN_LIBRARY_SRC} ${MASSV_LIBRARY_SRC})
add_dependencies(checkdmpart gem phy modelutils ${VGRID_LIBRARY_SRC} ${MPIEXTRAS_LIBRARY_SRC} ${RPNCOMM_LIBRARY_SRC} ${RMN_LIBRARY_SRC} ${MASSV_LIBRARY_SRC})

install(TARGETS maingemdm checkdmpart gemgrid prgemnml RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
