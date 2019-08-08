add_executable(mainyy2global mainyy2global.F90)
target_link_libraries(mainyy2global gem phy modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(mainyy2global gem phy modelutils vgrid rmn massv)

add_executable(flipit mainflipit.F90)
target_link_libraries(flipit modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(flipit modelutils vgrid rmn massv)

add_executable(r.fstinfo fstinfo.f)
add_executable(voir voir.f90)

set_target_properties(r.fstinfo PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(voir PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(r.fstinfo rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
target_link_libraries(voir rmn
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )

install(TARGETS mainyy2global flipit r.fstinfo voir RUNTIME
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/orig-scripts/scripts
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}
	USE_SOURCE_PERMISSIONS)
install(FILES ${CMAKE_SOURCE_DIR}/orig-scripts/runprep
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${CMAKE_SOURCE_DIR}/orig-scripts/runmod
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${ROOT}/orig-scripts/mybatchjob.pbs
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
install(FILES ${ROOT}/README.run
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION}
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
install(CODE "execute_process(COMMAND
        /bin/ln -sf ${ROOT}/../work/configurations ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION} ) ")
