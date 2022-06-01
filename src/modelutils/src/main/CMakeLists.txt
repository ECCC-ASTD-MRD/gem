message(STATUS "(EC) Generating modelutils main Makefile")

set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})

add_executable(yy2global mainyy2global.F90)
target_link_libraries(yy2global modelutils vgrid rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(yy2global modelutils vgrid rmn)

add_executable(flipit mainflipit.F90)
target_link_libraries(flipit modelutils vgrid rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(flipit modelutils vgrid rmn)

add_executable(yydecode mainyydecode.F90)
target_link_libraries(yydecode modelutils vgrid rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(yydecode modelutils vgrid rmn)

add_executable(yyencode mainyyencode.F90)
target_link_libraries(yyencode modelutils vgrid rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(yyencode modelutils vgrid rmn)

install(TARGETS yy2global flipit yydecode yyencode RUNTIME
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

