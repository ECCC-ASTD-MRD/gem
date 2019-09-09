add_definitions(-DAIX)

set(CMAKE_C_COMPILER "xlc")
set(CMAKE_Fortran_COMPILER "xlf_r")

set(MPI_C_COMPILER "xlc")
set(MPI_Fortran_COMPILER "xlf_r")

set(MPI_C_LIBRARIES "-L/usr/mpi/gcc/openmpi-4.0.0rc5/lib64 -lmpi")
set(MPI_C_INCLUDE_PATH "/usr/mpi/gcc/openmpi-4.0.0rc5/include")

set(MPI_Fortran_LIBRARIES "-L/usr/mpi/gcc/openmpi-4.0.0rc5/lib64 -lmpi")
set(MPI_Fortran_INCLUDE_PATH "/usr/mpi/gcc/openmpi-4.0.0rc5/include")

set(CMAKE_C_FLAGS "-Wl,--allow-shlib-undefined -Wtrigraphs -I. -qfunctrace  -q64 -qtbtable=full -qflttrap=ov:zerodivide:enable:imp -D_REENTRANT -D_THREAD_SAFE -qsmp=noauto -DWITH_OpenMP" CACHE STRING "C compiler flags" FORCE)

set(MPI_C_COMPILE_FLAGS "-Wl,--allow-shlib-undefined -Wtrigraphs -I. -qfunctrace -qcpluscmt -qweaksymbol -qlargepage -q64 -qtbtable=full -qflttrap=ov:zerodivide:enable:imp -mpi -D_REENTRANT -D_THREAD_SAFE -qsmp=noauto -DWITH_OpenMP" CACHE STRING "C compiler flags" FORCE)

set(CMAKE_Fortran_FLAGS "-qxlf2003=polymorphic -qnotrigraph -qstrict -qtbtable=full -q64 -qcache=auto -qarch=auto -qextname -qlargepage -qtune=auto -qnosave -qflttrap=ov:zerodivide:enable:imp -qsigtrap=xl__trce -qfloat=nofold -I. -O2 -qsmp=noauto" CACHE STRING "Fortran compiler flags" FORCE)
set(MPI_Fortran_COMPILE_FLAGS "-mpi" FORCE)

set(MPI_Fortran_COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}" CACHE STRING "Fortran compiler flags")

set(CMAKE_EXE_LINKER_FLAGS "-assume byterecl")

