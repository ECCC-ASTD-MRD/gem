message(STATUS "Starting generation of rpnphy Makefile")

add_definitions(-DECCCGEM)

include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})

file(GLOB FTNSRC */*.F90)

add_compile_options( "${EXTRA_FLAGS}" )

add_library (phy STATIC ${FTNSRC})

target_link_libraries(phy
   modelutils
)

target_include_directories(phy PRIVATE
   base
   include
   series
   surface
   utils
   ${ROOT}/src/modelutils/include
)
