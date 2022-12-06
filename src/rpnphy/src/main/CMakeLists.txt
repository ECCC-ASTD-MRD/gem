message(STATUS "(EC) Generating rpnphy binaries Makefile")

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   rpnphy
   modelutils
   vgrid
   rpncomm
   rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)
add_dependencies(feseri rpnphy ${modelutils_SRC} ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC})

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
   tdpack
   modelutils
   vgrid
   rpncomm
   rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
)
add_dependencies(prphynml rpnphy ${modelutils_SRC} ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC})

install(TARGETS feseri prphynml RUNTIME DESTINATION bin)
