message(STATUS "(EC) Generating gemdyn libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
   set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
endif()

file(GLOB F_FILES */*.F90)
file(GLOB C_FILES  */*.c)
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/gemdyn/src/main/main*" )
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/gemdyn/src/no_mpi/main*" )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules )

add_library(gemdyn STATIC ${F_FILES} ${C_FILES})
add_dependencies(gemdyn rpnphy modelutils vgrid tdpack mpi-extras rpncomm rmn)
target_include_directories(gemdyn PRIVATE
   ${CMAKE_CURRENT_SOURCE_DIR}/../include
   ${CMAKE_CURRENT_SOURCE_DIR}/base
   ${CMAKE_SOURCE_DIR}/src/modelutils/include
   ${CMAKE_BINARY_DIR}/src/tdpack
   ${CMAKE_BINARY_DIR}/src/vgrid
   ${CMAKE_BINARY_DIR}/src/modelutils/modelutils/modules
   ${CMAKE_BINARY_DIR}/src/rpnphy/rpnphy/modules
)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  if (EXTRA_CHECKS AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_compile_options(gemdyn PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:Fortran>>:-fcheck=bounds>" )
    message(STATUS "(EC) Extra debug flags for gemdyn: -fcheck=bounds")
  endif()
endif()

