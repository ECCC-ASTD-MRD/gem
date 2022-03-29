message(STATUS "(EC) Generating rpnphy libraries Makefile")

include(FindMPI OPTIONAL)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB F_FILES main/*.F90 api/*.F90 base/*.F90 series/*.F90 surface/*.F90 utils/*.F90)
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*" )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )

add_library(rpnphy STATIC ${F_FILES})
target_include_directories(rpnphy BEFORE PRIVATE 
  series
  surface )

include_directories(
  ${CMAKE_SOURCE_DIR}/src/modelutils/include
  ${CMAKE_BINARY_DIR}/src/vgrid
  ${CMAKE_BINARY_DIR}/src/modelutils/modelutils/modules )

target_include_directories(rpnphy PRIVATE ../include base utils)
add_dependencies(rpnphy rpnphy_cpl_stubs)
add_dependencies(rpnphy rpnphy_chm_stubs)

file(GLOB F_CHMSTUBS_FILES chm_stubs/*.F90)
add_library(rpnphy_chm_stubs STATIC ${F_CHMSTUBS_FILES})
add_library(rpnphy_cpl_stubs STATIC cpl_stubs/cpl_itf_stubs.F90)
