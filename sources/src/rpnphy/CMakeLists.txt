message(STATUS "Starting generation of rpnphy Makefile")

project(RPNPHY Fortran)

file(GLOB FTNSRC */*.F90)

include_directories(series)
include_directories(${ROOT}/src/modelutils/include)
include_directories(AFTER include base utils ${ROOT}/src/modelutils/include)

add_library (phy STATIC ${FTNSRC})
