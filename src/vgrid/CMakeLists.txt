message(STATUS "Starting generation of vgrid Makefile")

file(GLOB FTN_SOURCES *.F90)
file(GLOB CSRC  *.c)

message(STATUS "Sources C: ${CSRC}")
message(STATUS "Sources Fortran: ${FTN_SOURCES}")

add_library (vgrid STATIC ${FTN_SOURCES} ${CSRC})
