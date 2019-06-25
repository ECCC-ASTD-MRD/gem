message(STATUS "Starting generation of modelutils Makefile")
project(MODELUTILS Fortran)

include(FindMPI OPTIONAL)
#include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB CSRC  */*.c)
file(GLOB FTNSRC  */*.F90)

include_directories(AFTER base include)

add_compile_options( "${EXTRA_FLAGS}" )

add_library (modelutils STATIC ${FTNSRC} ${CSRC})
