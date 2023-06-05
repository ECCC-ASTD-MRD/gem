message(STATUS "(EC) Generating filetype Makefile")

add_executable(r.filetype filetyp.f90)
target_link_libraries(r.filetype PRIVATE rmn::rmn)

install(TARGETS r.filetype RUNTIME DESTINATION bin)
