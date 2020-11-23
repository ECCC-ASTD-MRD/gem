# CMAKE_<LANG>_COMPILER variables MUST contain the full path of the compiler for <LANG>
find_program(CMAKE_C_COMPILER icc)
find_program(CMAKE_Fortran_COMPILER ifort)

add_definitions(-DLittle_Endian)

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2")
set(CMAKE_C_FLAGS "-mkl -Wl,--allow-shlib-undefined -Wtrigraphs -xCORE-AVX512 -fpic -traceback -fp-model precise" CACHE STRING "C compiler flags" FORCE)

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -ftrapuv")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O2")
set(CMAKE_Fortran_FLAGS "-mkl -assume byterecl -convert big_endian -xCORE-AVX512 -align array32byte -fpe0 -fpic -traceback -ip -diag-disable=cpu-dispatch -static-intel -diag-disable 7713 -diag-disable 10212 -diag-disable 5140 -fp-model source" CACHE STRING "Fortran compiler flags" FORCE)

set(CMAKE_EXE_LINKER_FLAGS_INIT "-Wl,--allow-shlib-undefined,--allow-multiple-definition -fpic -static-intel -mkl")
#set(CMAKE_EXE_LINKER_FLAGS_INIT "-fpic -static-intel -mkl")
