message(STATUS "Starting generation of gemdyn Makefile")
project(GEMDYN Fortran)

include(FindMPI OPTIONAL)
include_directories(${MPI_INCLUDE_PATH})
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC */*.F90)
file(GLOB CSRC  */*.c)

include_directories(${ROOT}/../src/rpnphy/include
                    ${ROOT}/../src/modelutils/include)

include_directories(AFTER include base)

if ( "${CMAKE_C_COMPILER}" MATCHES ".*/icc" )
  add_compile_options( "-std08" )
  message( STATUS "Adding -std08 when compiler = " ${CMAKE_C_COMPILER} )
elseif( "${CMAKE_C_COMPILER}" MATCHES ".*/gcc" )
  add_compile_options( "-std=f2008" )
  add_compile_options( "-pedantic" )
  message( STATUS "Adding -std=f2008 -pedantic when compiler = " ${CMAKE_C_COMPILER} )
endif()

add_library (gem STATIC ${FTNSRC} ${CSRC})
