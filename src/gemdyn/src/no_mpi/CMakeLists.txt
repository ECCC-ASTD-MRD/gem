message(STATUS "(EC) Generating gem_monitor* and toc2nml makefiles")

add_executable(toc2nml maintoc2nml.F90)
target_link_libraries(toc2nml
   gemdyn
   rpnphy
   modelutils
   vgrid
#   tdpack
   rmn
   ${ompfftw}
   fftw3
   ${LAPACK_LIBRARIES}
   ${BLAS_LIBRARIES}
   ${EXTRA_LIBRARIES}
   )
target_include_directories(toc2nml PRIVATE
  ${CMAKE_BINARY_DIR}/src/vgrid
 )

add_dependencies(toc2nml gemdyn rpnphy modelutils vgrid rmn)

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
set_target_properties(gem_monitor_output PROPERTIES LINKER_LANGUAGE C)
set_target_properties(gem_monitor_end PROPERTIES LINKER_LANGUAGE C)

install(TARGETS toc2nml gem_monitor_output gem_monitor_end RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

