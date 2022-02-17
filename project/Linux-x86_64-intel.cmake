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

set(CMAKE_C_FLAGS "-Wl,--allow-shlib-undefined -Wtrigraphs -fpic -traceback -fp-model precise" CACHE STRING "C compiler flags" FORCE)
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O2")

set(CMAKE_Fortran_FLAGS "-assume byterecl -convert big_endian -fpe0 -fpic -traceback -static-intel -diag-disable 7713 -diag-disable 10212 -diag-disable 5268 -diag-disable 5140 -fp-model source" CACHE STRING "Fortran compiler flags" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "-g")
set(CMAKE_Fortran_FLAGS_RELEASE "-O2")

set(CMAKE_EXE_LINKER_FLAGS_INIT "-fpic -static-intel")

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2021)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qmkl")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl")
  set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -qmkl")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
  set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -mkl")
endif()

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

execute_process(COMMAND sh "-c" "${CMAKE_CURRENT_SOURCE_DIR}/intel_find_best_arch.sh" OUTPUT_VARIABLE best_arch ERROR_QUIET)
message(STATUS "Most performant -march for the Intel compilers on this hardware: ${best_arch}")
message(STATUS "You could add this flag to CMAKE_C_FLAGS")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=${best_arch}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -march=${best_arch}")
