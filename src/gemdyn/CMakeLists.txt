message(STATUS "Starting generation of gemdyn Makefile")

include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB FTNSRC */*.F90)
file(GLOB CSRC  */*.c)
list(FILTER FTNSRC EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/gemdyn/main/main*" )

add_library (gem STATIC ${FTNSRC} ${CSRC})
target_include_directories(gem PRIVATE
   ${CMAKE_SOURCE_DIR}/src/rpnphy/include
   ${CMAKE_SOURCE_DIR}/src/modelutils/include
   ${CMAKE_BINARY_DIR}/src/rpnphy
   ${CMAKE_BINARY_DIR}/src/modelutils
   ${CMAKE_BINARY_DIR}/src/vgrid
   include
   base
)

# Adding specific flags for GEM
include(${CMAKE_SOURCE_DIR}/src/flags-GEM.cmake)

if (NOT "${CMAKE_SYSTEM_NAME}" STREQUAL "CrayLinuxEnvironment")
  set(EXTRA_LIBRARIES "-lhwloc")
endif()

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  if (EXTRA_DEBUG AND "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    target_compile_options( gem PRIVATE "$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:Fortran>>:-fcheck=bounds>" )
    message(STATUS "Extra debug flags for gemdyn: -fcheck=bounds")
  endif()
endif()

add_subdirectory(main gemdyn-main)
add_subdirectory(no_mpi gemdyn_no_mpi)
