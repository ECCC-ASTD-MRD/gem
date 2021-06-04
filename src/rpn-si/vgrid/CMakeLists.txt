message(STATUS "Starting generation of vgrid Makefile")

file(GLOB FTN_SOURCES *.F90)
file(GLOB CSRC  *.c)

add_library (vgrid STATIC ${FTN_SOURCES} ${CSRC})
