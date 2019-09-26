message(STATUS "Starting generation of fststat Makefile")

add_executable(fststat fststatm.f fststat.f statfld4.f)
target_link_libraries(fststat rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_dependencies(fststat rmn)

install(TARGETS fststat RUNTIME DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
