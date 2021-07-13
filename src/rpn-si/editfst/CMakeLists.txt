message(STATUS "Starting generation of editfst Makefile")

FILE(GLOB FTN_SOURCES *.f90)

add_executable(editfst bidon.f90 ${FTN_SOURCES})
target_link_libraries(editfst
   rmn
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
)

add_dependencies(editfst rmn)

install(TARGETS editfst RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
