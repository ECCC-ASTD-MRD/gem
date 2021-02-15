message(STATUS "Starting generation of gemdyn Makefile")

include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB FTNSRC */*.F90)
file(GLOB CSRC  */*.c)

add_compile_options( "${EXTRA_FLAGS}" )

add_library (gem STATIC ${FTNSRC} ${CSRC})
target_include_directories(gem PRIVATE
   ${ROOT}/src/rpnphy/include
   ${ROOT}/src/modelutils/include
   ${CMAKE_BINARY_DIR}/src/rpnphy
   ${CMAKE_BINARY_DIR}/src/modelutils
   ${CMAKE_BINARY_DIR}/src/vgrid
   include
   base
)

