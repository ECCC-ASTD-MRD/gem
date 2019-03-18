add_executable(maingemdm.Abs maingemdm.F90)
target_link_libraries(maingemdm.Abs gem phy modelutils vgrid rpncomm rmn massv
                      ${MPI_Fortran_LIBRARIES} 
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(maingemdm.Abs gem phy modelutils vgrid rpncomm rmn massv)

add_executable(maingemgrid.Abs maingemgrid.F90)
target_link_libraries(maingemgrid.Abs gem phy modelutils vgrid rpncomm rmn massv
                      ${MPI_Fortran_LIBRARIES} 
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(maingemgrid.Abs gem phy modelutils vgrid rpncomm rmn massv)

add_executable(mainyy2global.Abs mainyy2global.F90)
target_link_libraries(mainyy2global.Abs gem phy modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(mainyy2global.Abs gem phy modelutils vgrid rmn massv)

add_executable(flipit.Abs mainflipit.F90)
target_link_libraries(flipit.Abs modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(flipit.Abs modelutils vgrid rmn massv)

add_executable(r.fstinfo fstinfo.f)
add_executable(voir voir.f90)

set_target_properties(r.fstinfo PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(voir PROPERTIES LINKER_LANGUAGE Fortran)

target_link_libraries(r.fstinfo rmn 
                      ${Fortran_LIBRARIES}
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
target_link_libraries(voir rmn 
                      ${Fortran_LIBRARIES}
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )

install(TARGETS maingemdm.Abs maingemgrid.Abs mainyy2global.Abs flipit.Abs r.fstinfo voir RUNTIME 
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
install(FILES ${ROOT}/README.run
        DESTINATION ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION} 
        PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
install(CODE "execute_process(COMMAND
        ${CMAKE_SOURCE_DIR}/lncfg.sh ${ROOT}/../work/configurations ${ROOT}/../work/work-${OS}-${COMPILER}-${COMP_VERSION} ) ")
