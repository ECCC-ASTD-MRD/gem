message(STATUS "Starting generation of fstinfo Makefile")

add_executable(r.fstinfo fstinfo.f)
set_target_properties(r.fstinfo PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(r.fstinfo rmn )

install(TARGETS r.fstinfo  RUNTIME
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

