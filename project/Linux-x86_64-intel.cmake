# CMAKE_<LANG>_COMPILER variables MUST contain the full path of the compiler for <LANG>
find_program(CMAKE_C_COMPILER icc)
find_program(CMAKE_Fortran_COMPILER ifort)

# I don't know why, but enable_language exmpties CMAKE_BUILD_TYPE!
# We therefore have to back it up and restore it after enable_language
set(TMP_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Enable the two languages that are used
enable_language(C)
enable_language(Fortran)

# Reset CMAKE_BUILD_TYPE
set(CMAKE_BUILD_TYPE ${TMP_BUILD_TYPE})

# find_package() commands can only be called after the languages have been 
# eneabled or they will fail

add_definitions(-DLittle_Endian -DWITH_intel)

set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")
set(CMAKE_C_FLAGS "-mkl -Wl,--allow-shlib-undefined -Wtrigraphs -fpic -traceback -fp-model precise ${compFlags}" CACHE STRING "C compiler flags" FORCE)

set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
set(CMAKE_Fortran_FLAGS "-mkl -assume byterecl -convert big_endian -fpe0 -fpic -reentrancy threaded -traceback -threads -static-intel -diag-disable 7713 -diag-disable 10212 -diag-disable 5140 -fp-model source ${compFlags}" CACHE STRING "Fortran compiler flags" FORCE)

set(CMAKE_EXE_LINKER_FLAGS_INIT "-fpic -static-intel -mkl")

# The Intel compilers do not support -march=native like GCC does
# This will have to be edited in order to get the most out of the hardware

# If you run into unknown instruction errors at runtime, you will have to change the
# target architecture specified for "-march".

# Here are architecture names given in descending order of performance:
#     skylake-avx512
#     skylake
#     core-avx2
#     core-avx-i
#     corei7-avx
#     core2

# Set flags according to the compiler version
# We assume Intel 19 is used on skylake-avx512 or newer hardware
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "19")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=skylake-avx512")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=skylake-avx512")

elseif(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "16")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=core-avx-i")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=core-avx-i")

else()
   message(FATAL_ERROR "Unknown Intel compiler version!  You can probably make it work by editing ${CMAKE_CURRENT_LIST_FILE}")
endif()



