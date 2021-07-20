# CMAKE_<LANG>_COMPILER variables MUST contain the full path of the compiler for <LANG>
find_program(CMAKE_C_COMPILER icc)
find_program(CMAKE_Fortran_COMPILER ifort)

# I don't know why, but enable_language empties CMAKE_BUILD_TYPE!
# We therefore have to back it up and restore it after enable_language
set(TMP_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Enable the two languages that are used
enable_language(C)
enable_language(Fortran)

# Reset CMAKE_BUILD_TYPE
set(CMAKE_BUILD_TYPE ${TMP_BUILD_TYPE})

# find_package() commands can only be called after the languages have been 
# enabled or they will fail

add_definitions(-DLittle_Endian)

set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS "-mkl -Wtrigraphs -fpic -traceback -fp-model precise" CACHE STRING "C compiler flags" FORCE)

set(CMAKE_Fortran_FLAGS_DEBUG "-g -ftrapuv")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS "-mkl -assume byterecl -convert big_endian -fpe0 -fpic -traceback -static-intel -diag-disable 7713 -diag-disable 10212 -diag-disable 5140 -fp-model source" CACHE STRING "Fortran compiler flags" FORCE)

set(CMAKE_EXE_LINKER_FLAGS_INIT "-fpic -static-intel -mkl")

set(EXTRA_LIBRARIES "-lopen-pal")
