message(STATUS "Starting generation of fstcomp Makefile")

file(GLOB FTNSRC LIST_DIRECTORIES false RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.F90 *.F *.f *.f90)

add_executable(fstcomp ${FTNSRC})
target_link_libraries(fstcomp
   rmn
)

add_dependencies(fstcomp rmn)

install(TARGETS fstcomp RUNTIME DESTINATION ${ROOT}/work-${OS}-${COMPILER}-${COMP_VERSION}/bin)
