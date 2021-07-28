message(STATUS "Starting generation of rpnphy Makefile")

project(RPNPHY Fortran)

include(FindMPI OPTIONAL)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC */*.F90)
list(FILTER FTNSRC EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*" )

include_directories(series)
include_directories(surface)
include_directories(${CMAKE_SOURCE_DIR}/src/modelutils/include)
include_directories(${CMAKE_BINARY_DIR}/src/modelutils)
include_directories(${CMAKE_BINARY_DIR}/src/vgrid)
include_directories(AFTER include base utils ${CMAKE_SOURCE_DIR}/src/modelutils/include)

add_library (phy STATIC ${FTNSRC})

add_subdirectory(main rpnphy-main)
