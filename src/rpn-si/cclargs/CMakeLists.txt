message(STATUS "Starting generation of cclargs_lite Makefile")

add_definitions(-DNOUI)

add_executable(cclargs_lite cclargs_lite.c)

install(TARGETS cclargs_lite RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
