message(STATUS "(EC) Generating modelutils Makefile")

enable_language(C Fortran)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/modelutils_version.in ${CMAKE_CURRENT_SOURCE_DIR}/include/modelutils_version.inc @ONLY)

add_subdirectory(src/main modelutils-main)
add_subdirectory(src modelutils)
