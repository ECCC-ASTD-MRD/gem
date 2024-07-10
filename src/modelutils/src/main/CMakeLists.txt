message(STATUS "(EC) Generating modelutils main Makefile")

add_executable(yy2global yy2global.F90)
set_target_properties(yy2global PROPERTIES COMPILE_DEFINITIONS YY2GLOBAL_VERSION="${YY2GLOBAL_VERSION}")
target_link_libraries(yy2global
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
 
add_executable(flipit flipit.F90)
set_target_properties(flipit PROPERTIES COMPILE_DEFINITIONS FLIPIT_VERSION="${FLIPIT_VERSION}")
target_link_libraries(flipit
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  )

add_executable(yydecode yydecode.F90)
set_target_properties(yydecode PROPERTIES COMPILE_DEFINITIONS YYDECODE_VERSION="${YYDECODE_VERSION}")
target_link_libraries(yydecode
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} )

add_executable(yyencode yyencode.F90)
set_target_properties(yyencode PROPERTIES COMPILE_DEFINITIONS YYENCODE_VERSION="${YYENCODE_VERSION}")
target_link_libraries(yyencode
  modelutils
  vgrid::vgrid
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} )

add_executable(test_integrals maintest_integrals.F90)
target_link_libraries(test_integrals
  modelutils
  rmn::rmn
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES} )

install(TARGETS yy2global flipit yydecode yyencode test_integrals RUNTIME DESTINATION bin)
