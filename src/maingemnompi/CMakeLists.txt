message(STATUS "Starting generation of gem_monitor* Makefile")

set(CMAKE_Fortran_COMPILER ${Fortran_COMPILER})
set(CMAKE_C_COMPILER ${C_COMPILER})

if (WITH_OPENMP)
   set(ompfftw "-lfftw3_omp")
endif()

add_executable(gem_monitor_output gem_monitor_output.c)
add_executable(gem_monitor_end gem_monitor_end.c)
set_target_properties(gem_monitor_output PROPERTIES LINKER_LANGUAGE C)
set_target_properties(gem_monitor_end PROPERTIES LINKER_LANGUAGE C)

install(TARGETS gem_monitor_output gem_monitor_end RUNTIME 
        DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)

