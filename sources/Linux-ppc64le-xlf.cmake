add_definitions(-DLittle_Endian -DWITH_OpenMP -D_REENTRANT -D_THREAD_SAFE)

set(CMAKE_C_COMPILER "xlc_r")
set(CMAKE_Fortran_COMPILER "xlf_r")

set(MPI_C_COMPILER "mpixlc")
set(MPI_Fortran_COMPILER "mpixlf")

set(CMAKE_C_FLAGS "-Wl,--allow-shlib-undefined -Wtrigraphs -qflttrap=ov:zerodivide:enable:imp -O2 -qtune=auto -qsmp=omp" CACHE STRING "C compiler flags" FORCE)

set(MPI_C_COMPILE_FLAGS "-Wl,--allow-shlib-undefined -Wtrigraphs -qflttrap=ov:zerodivide:enable:imp -qsmp=omp" CACHE STRING "C compiler flags" FORCE)

set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic -qstrict -qcache=auto -qarch=auto -qextname -qtune=auto -qnosave -qflttrap=ov:zerodivide:enable:imp -qfloat=nofold -O2 -qsmp=omp" CACHE STRING "Fortran compiler flags" FORCE)
set(MPI_Fortran_COMPILE_FLAGS "-mpi" FORCE)

set(MPI_Fortran_COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}" CACHE STRING "Fortran compiler flags")

set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-qmkshrobj")
