message(STATUS "(EC) Generating rpnphy Makefile")

project(RPNPHY Fortran)

add_subdirectory(main rpnphy-main)

include(FindMPI OPTIONAL)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC */*.F90)
list(FILTER FTNSRC EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*" )

add_library (phy STATIC ${FTNSRC})
target_include_directories(phy BEFORE PRIVATE 
  series
  surface
  ${CMAKE_SOURCE_DIR}/src/modelutils/include
  ${CMAKE_BINARY_DIR}/src/modelutils/modelutils/modules
  ${VGRID_INCLUDE_DIR}
  )
target_include_directories(phy PRIVATE include base utils)

