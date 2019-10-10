message(STATUS "Starting generation of fstcomp Makefile")

add_executable(fstcomp statfldx.f90 fstcomp.f)
target_link_libraries(fstcomp
   rmn
)

add_dependencies(fstcomp rmn)

install(TARGETS fstcomp RUNTIME DESTINATION ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
