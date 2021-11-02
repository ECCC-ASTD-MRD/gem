message(STATUS "(EC) Generating filetype Makefile")

add_executable(r.filetype filetype.f90)
set_target_properties(r.filetype PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(r.filetype rmn )

install(TARGETS r.filetype RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

