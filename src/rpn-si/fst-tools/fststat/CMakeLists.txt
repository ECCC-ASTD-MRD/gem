message(STATUS "(EC) Generating fststat Makefile")

add_executable(fststat fststatm.f fststat.f statfld4.f)
target_link_libraries(fststat rmn)

add_dependencies(fststat rmn)

install(TARGETS fststat RUNTIME DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
