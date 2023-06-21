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
add_dependencies(feseri rpnphy ${modelutils_SRC})

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
   tdpack::tdpack 
   modelutils
   vgrid::vgrid
   rpncomm::rpncomm
   rmn::rmn-ompi
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)
add_dependencies(prphynml rpnphy ${modelutils_SRC})

install(TARGETS feseri prphynml RUNTIME DESTINATION bin)
