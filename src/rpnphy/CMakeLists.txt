message(STATUS "(EC) Generating rpnphy Makefile")

enable_language(C Fortran)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/rpnphy_version.in ${CMAKE_CURRENT_SOURCE_DIR}/include/rpnphy_version.inc @ONLY)

add_subdirectory(src/main rpnphy-main)
add_subdirectory(src rpnphy)

