add_definitions(-D__GNUC__=4 -DAMD64 -DLINUX_X86_64 -DLittle_Endian)

set(CMAKE_C_FLAGS "-g -D_REENTRANT" CACHE STRING "C compiler flags" FORCE)
set(CMAKE_Fortran_FLAGS "-g -h byteswapio -em -J ${CMAKE_CURRENT_BINARY_DIR}/${BUILD}/modules" CACHE STRING "Fortran compiler flags" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-s -fopenmp -fpic")
set(MPI_Fortran_COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}" CACHE STRING "Fortran compiler flags")

set(LAPACK_LIBRARIES "-llapack")
set(BLAS_LIBRARIES "-lblas")
