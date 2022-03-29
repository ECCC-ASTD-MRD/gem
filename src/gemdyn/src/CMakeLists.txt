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

target_include_directories(gemdyn PRIVATE
   ${CMAKE_CURRENT_SOURCE_DIR}/../include
   ${CMAKE_CURRENT_SOURCE_DIR}/base
   ${CMAKE_SOURCE_DIR}/src/modelutils/include
   ${CMAKE_BINARY_DIR}/src/vgrid
   ${CMAKE_BINARY_DIR}/src/modelutils/modelutils/modules
   ${CMAKE_BINARY_DIR}/src/rpnphy/rpnphy/modules
)

