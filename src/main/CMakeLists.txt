add_executable(yy2global mainyy2global.F90)
target_link_libraries(yy2global gem phy modelutils vgrid rmn massv
                      ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${EXTRA_LIBRARIES} )
add_dependencies(yy2global gem phy modelutils vgrid rmn massv)

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

install(TARGETS yy2global flipit r.fstinfo voir RUNTIME
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)


