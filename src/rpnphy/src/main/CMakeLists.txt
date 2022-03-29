message(STATUS "(EC) Generating rpnphy binaries Makefile")

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   rpnphy
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
add_dependencies(feseri rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
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
add_dependencies(prphynml rpnphy modelutils vgrid mpi-extras rpncomm rmn massv)

install(TARGETS feseri prphynml RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

