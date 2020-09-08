add_executable(mainyy2global mainyy2global.F90)
target_link_libraries(mainyy2global gem phy modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(mainyy2global gem phy modelutils vgrid rmn massv)

add_executable(flipit mainflipit.F90)
target_link_libraries(flipit modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(flipit modelutils vgrid rmn massv)

#add_executable(r.fstinfo fstinfo.f)
#add_executable(voir voir.f90)

#set_target_properties(r.fstinfo PROPERTIES LINKER_LANGUAGE Fortran)
#set_target_properties(voir PROPERTIES LINKER_LANGUAGE Fortran)

#target_link_libraries(r.fstinfo rmn
#                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
#target_link_libraries(voir rmn
#                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )

#install(TARGETS mainyy2global flipit r.fstinfo voir RUNTIME
#        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

install(TARGETS mainyy2global flipit  RUNTIME
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
