message(STATUS "Starting generation of rpnphy Makefile")

project(RPNPHY Fortran)

include(FindMPI OPTIONAL)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC */*.F90)

# Adding specific flags for GEM
include(${CMAKE_SOURCE_DIR}/src/flags-GEM.cmake)

include_directories(series)
include_directories(${CMAKE_SOURCE_DIR}/src/modelutils/include)
include_directories(${CMAKE_BINARY_DIR}/src/modelutils)
include_directories(${CMAKE_BINARY_DIR}/src/vgrid)
include_directories(AFTER include base utils ${CMAKE_SOURCE_DIR}/src/modelutils/include)

add_library (phy STATIC ${FTNSRC})
