message(STATUS "(EC) Generating rpnphy binaries Makefile")

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   rpnphy
   modelutils
   vgrid::vgrid
   rpncomm::rpncomm
   rmn::rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
   tdpack
   modelutils
   vgrid::vgrid
   rpncomm::rpncomm
   rmn::rmn-ompi
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)

install(TARGETS feseri prphynml RUNTIME DESTINATION bin)
