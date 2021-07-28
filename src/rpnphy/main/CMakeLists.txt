message(STATUS "Starting generation of rpnphy main Makefile")

set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(feseri phy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(prphynml phy modelutils vgrid mpi-extras rpncomm rmn massv)

install(TARGETS feseri prphynml RUNTIME 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
