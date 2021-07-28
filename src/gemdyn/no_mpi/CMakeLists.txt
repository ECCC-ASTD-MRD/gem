message(STATUS "Starting generation of gem_monitor Makefile")

set(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})

add_executable(toc2nml toc2nml.F90)
target_link_libraries(toc2nml
   gem
   phy
   modelutils
   vgrid
   mpi-extras
   rpncomm
   rmn
   massv
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
   )
target_include_directories(toc2nml PRIVATE
   ${CMAKE_BINARY_DIR}/src/vgrid
)
add_dependencies(toc2nml gem phy modelutils vgrid mpi-extras rpncomm rmn massv)

add_executable(gem_monitor_output maingem_monitor_output.c)
add_executable(gem_monitor_end maingem_monitor_end.c)
set_target_properties(gem_monitor_output PROPERTIES LINKER_LANGUAGE C)
set_target_properties(gem_monitor_end PROPERTIES LINKER_LANGUAGE C)

install(TARGETS toc2nml gem_monitor_output gem_monitor_end RUNTIME 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

