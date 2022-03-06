message(STATUS "Starting generation of rpncomm Makefile")

add_definitions (-DMARKER=RPN_COMM_marker_4542)

find_package(MPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB CSRC *.c)
file(GLOB FTNSRC *.F90)

add_library (rpncomm STATIC ${FTNSRC} ${CSRC})

target_compile_options( rpncomm PRIVATE $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-w> )
