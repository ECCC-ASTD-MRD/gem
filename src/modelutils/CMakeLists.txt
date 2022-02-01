message(STATUS "Starting generation of modelutils Makefile")

include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB CSRC  */*.c)
file(GLOB FTNSRC  */*.F90)

add_library (modelutils STATIC ${FTNSRC} ${CSRC})
target_include_directories(modelutils PRIVATE base include)
include_directories( ${CMAKE_BINARY_DIR}/src/vgrid )
