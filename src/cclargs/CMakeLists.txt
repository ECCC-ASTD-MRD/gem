message(STATUS "Starting generation of cclargs_lite Makefile")

project(CCLARGS_LITE C)

add_definitions(-DNOUI)

add_executable(cclargs_lite cclargs_lite.c)

#target_link_libraries(cclargs_lite
#                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )

#set_target_properties(cclargs_lite PROPERTIES LINKER_LANGUAGE C)

install(TARGETS cclargs_lite RUNTIME DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
