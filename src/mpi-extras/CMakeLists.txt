message(STATUS "Starting generation of mpi-extras Makefile")

find_package(MPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC RPN_MPI*.F90)

add_library (mpi-extras STATIC ${FTNSRC})
