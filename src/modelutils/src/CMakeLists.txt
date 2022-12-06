message(STATUS "(EC) Generating modelutils librairies Makefile")

find_package(MPI)
if (MPI_FOUND)
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

file(GLOB C_FILES */*.c)
file(GLOB F_FILES main/flipit.F90 main/yy2global.F90 main/yydecode.F90 main/yyencode.F90 base/*F90 driver/*.F90 utils/*.F90)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )

add_library(modelutils STATIC ${F_FILES} ${C_FILES})
add_dependencies(modelutils vgrid tdpack)

target_include_directories(modelutils PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/base
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CMAKE_BINARY_DIR}/src/tdpack
  ${CMAKE_BINARY_DIR}/src/vgrid
)

add_library(modelutils_ov_ifort STATIC ov_ifort/ifort_fpe_handler.F90)
add_dependencies(modelutils_ov_ifort modelutils)

add_library(modelutils_stubs STATIC stubs/drv_dyn_itf_stubs.F90)
target_include_directories(modelutils_stubs PRIVATE
  ${CMAKE_BINARY_DIR}/src/vgrid
)
add_dependencies(modelutils_stubs modelutils)

add_library(modelutils_tmg_stubs STATIC tmg_stubs/tmg_stubs.F90)
