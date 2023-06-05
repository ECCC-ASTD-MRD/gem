message(STATUS "(EC) Generating gem_monitor* and toc2nml makefiles")

add_executable(toc2nml maintoc2nml.F90)
target_link_libraries(toc2nml PUBLIC
   gemdyn
   ${phy_LIBRARIES}
   modelutils
   vgrid::vgrid
   rmn::rmn
   ${FFTW}
   ${FFTW_DOUBLE_LIB}
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   )

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/maingem_monitor_output.c
  COMMAND sed 's/main_gem_monitor_output/main/' ${CMAKE_CURRENT_SOURCE_DIR}/gem_monitor_output.c > ${CMAKE_CURRENT_SOURCE_DIR}/maingem_monitor_output.c
  COMMENT "Generating maingem_monitor_output.c"
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gem_monitor_output.c
  USES_TERMINAL
  )

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/maingem_monitor_end.c
  COMMAND sed 's/main_gem_monitor_end/main/' ${CMAKE_CURRENT_SOURCE_DIR}/gem_monitor_end.c > ${CMAKE_CURRENT_SOURCE_DIR}/maingem_monitor_end.c
  COMMENT "Generating maingem_monitor_end.c"
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gem_monitor_end.c
  USES_TERMINAL
  )

add_executable(gem_monitor_output maingem_monitor_output.c)
add_executable(gem_monitor_end maingem_monitor_end.c)

install(TARGETS toc2nml gem_monitor_output gem_monitor_end RUNTIME DESTINATION bin)
