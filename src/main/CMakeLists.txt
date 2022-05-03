message(STATUS "(EC) Generating rpnphy main Makefile")

add_executable(feseri mainfeseri.F90)
target_link_libraries(feseri
   rpnphy
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(feseri rpnphy ${modelutils_SRC} ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC} ${massv_SRC})

add_executable(prphynml mainprphynml.F90)
target_link_libraries(prphynml
   rpnphy
   tdpack
   modelutils
   vgrid
   rpncomm
   rmn
   massv
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)
add_dependencies(prphynml rpnphy ${modelutils_SRC} ${vgrid_SRC} ${rpncomm_SRC} ${rmn_SRC} ${massv_SRC})

install(TARGETS feseri prphynml RUNTIME DESTINATION bin)
