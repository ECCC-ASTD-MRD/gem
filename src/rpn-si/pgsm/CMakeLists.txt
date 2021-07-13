message(STATUS "Starting generation of pgsm Makefile")

add_executable(pgsm pgsm.F90 f_pgsm.f90 c_pgsm.c)
target_link_libraries(pgsm
   rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_dependencies(pgsm rmn)

install(TARGETS pgsm RUNTIME DESTINATION DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
