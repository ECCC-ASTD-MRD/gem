find_program(CMAKE_C_COMPILER pgcc)
find_program(CMAKE_Fortran_COMPILER pgfortran)

# I don't know why, but enable_language empties CMAKE_BUILD_TYPE!
# We therefore have to back it up and restore it after enable_language
set(TMP_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Enable the two languages that are used
enable_language(C)
enable_language(Fortran)

# Reset CMAKE_BUILD_TYPE
set(CMAKE_BUILD_TYPE ${TMP_BUILD_TYPE})

# find_package() commands can only be called after the languages have been 
# eneabled or they will fail
add_definitions(-DLittle_Endian)

set(CMAKE_C_FLAGS "-g -fpic -I." CACHE STRING "C compiler flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(CMAKE_Fortran_FLAGS "-g -byteswapio -fast -Mvect=fuse,simd -Kieee -fpic -traceback" CACHE STRING "Fortran compiler flags" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")

#set(MPI_Fortran_COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS} ${CMAKE_Fortran_FLAGS}" CACHE STRING "MPI Fortran compiler flags")

set(CMAKE_EXE_LINKER_FLAGS_INIT "-Wl,--allow-shlib-undefined -fpic" CACHE STRING "Linker flags" FORCE)

set(LAPACK_LIBRARIES "-llapack")
set(BLAS_LIBRARIES "-lblas")
