message(STATUS "Starting generation of massv Makefile")

file(GLOB FTN_SOURCES *.f)

add_library (massv STATIC ${FTN_SOURCES})
