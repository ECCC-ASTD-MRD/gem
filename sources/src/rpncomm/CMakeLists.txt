message(STATUS "Starting generation of rpncomm Makefile")

project(RPNCOMM Fortran)

add_definitions (-DMARKER=RPN_COMM_marker_4051612)
include(FindMPI OPTIONAL)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB CSRC  *.c)
file(GLOB FTNSRC  *.F90)

add_library (rpncomm STATIC ${FTNSRC} ${CSRC})
