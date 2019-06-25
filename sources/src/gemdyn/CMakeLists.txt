message(STATUS "Starting generation of gemdyn Makefile")
project(GEMDYN Fortran)

include(FindMPI OPTIONAL)
#include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB FTNSRC */*.F90)
file(GLOB CSRC  */*.c)

include_directories(${ROOT}/src/rpnphy/include
                    ${ROOT}/src/modelutils/include)

include_directories(AFTER include base)
add_compile_options( "${EXTRA_FLAGS}" )

add_library (gem STATIC ${FTNSRC} ${CSRC})
