message(STATUS "(EC) Generating rpnphy main Makefile")

set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
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

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
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

add_dependencies(prphynml phy modelutils ${VGRID_LIBRARY_SRC} ${MPIEXTRAS_LIBRARY_SRC} ${RPNCOMM_LIBRARY_SRC} ${RMN_LIBRARY_SRC} ${MASSV_LIBRARY_SRC})
add_dependencies(feseri phy modelutils ${VGRID_LIBRARY_SRC} ${MPIEXTRAS_LIBRARY_SRC} ${RPNCOMM_LIBRARY_SRC} ${RMN_LIBRARY_SRC} ${MASSV_LIBRARY_SRC})

install(TARGETS feseri prphynml RUNTIME 
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
