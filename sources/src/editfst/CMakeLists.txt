message(STATUS "Starting generation of editfst Makefile")

project(EDITFST Fortran)

FILE(GLOB FTN_SOURCES *.f90)

add_executable(editfst bidon.f90 ${FTN_SOURCES})
target_link_libraries(editfst rmn
                      ${Fortran_LIBRARIES}
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )

add_dependencies(editfst rmn)

install(TARGETS editfst RUNTIME DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

