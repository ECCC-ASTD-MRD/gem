message(STATUS "(EC) Generating cclargs_lite Makefile")

add_definitions(-DNOUI)

add_executable(cclargs_lite cclargs_lite.c)

install(TARGETS cclargs_lite RUNTIME DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
