message(STATUS "(EC) Generating modelutils main Makefile")

add_executable(yy2global mainyy2global.F90)
target_link_libraries(yy2global
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  )

add_executable(flipit mainflipit.F90)
target_link_libraries(flipit
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  )

add_executable(yydecode mainyydecode.F90)
target_link_libraries(yydecode
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} )

add_executable(yyencode mainyyencode.F90)
target_link_libraries(yyencode
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} )

install(TARGETS yy2global flipit yydecode yyencode RUNTIME DESTINATION bin)
