find_program(CMAKE_C_COMPILER "icc")
find_program(CMAKE_Fortran_COMPILER "ifort")

find_program(MPI_C_COMPILER "cc")
find_program(MPI_Fortran_COMPILER "ftn")

# I don't know why, but enable_language empties CMAKE_BUILD_TYPE!
# We therefore have to back it up and restore it after enable_language
set(TMP_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Enable the two languages that are used
enable_language(C)
enable_language(Fortran)

# Reset CMAKE_BUILD_TYPE
set(CMAKE_BUILD_TYPE ${TMP_BUILD_TYPE})

add_definitions(-DLittle_Endian -DECCCGEM)

set(EXTRA_LIBRARIES "-liomp5 -L/opt/gcc/6.1.0/snos/lib64 -lcilkrts")

set(AVG_FLAG -O2)
set(FAST_FLAG -O2)
#set(FAST_FLAG -O3 -fast-transcendentals -no-prec-div -ip -no-prec-sqrt)

set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS "-traceback -fpic -xCORE-AVX512 -fp-model precise" CACHE STRING "C compiler flags" FORCE)

set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS "-traceback -assume byterecl -convert big_endian -xCORE-AVX512 -align array32byte -fpe0 -fpic -ip  -diag-disable=cpu-dispatch -static-intel -diag-disable 7713 -diag-disable 10212 -diag-disable 5140  -fp-model source -I." CACHE STRING "Fortran compiler flags" FORCE)

set(CMAKE_EXE_LINKER_FLAGS_INIT "-Wl,--allow-shlib-undefined,--allow-multiple-definition -fpic -Bstatic -static-libgcc -static")

