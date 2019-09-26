message(STATUS "Starting generation of rpnphy Makefile")

add_definitions(-DECCCGEM)

include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC */*.F90)

add_compile_options( "${EXTRA_FLAGS}" )

add_library (phy STATIC ${FTNSRC})
target_include_directories(phy PRIVATE
   series
   ${ROOT}/../src/modelutils/include
   include
   base
   utils
   ${ROOT}/../src/modelutils/include
)
