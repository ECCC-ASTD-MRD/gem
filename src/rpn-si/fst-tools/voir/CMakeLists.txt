message(STATUS "(EC) Generating voir Makefile")

add_executable(voir voir.f90)

set_target_properties(voir PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(voir rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )

install(TARGETS voir RUNTIME
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
