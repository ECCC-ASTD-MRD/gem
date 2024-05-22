cmake_minimum_required(VERSION 3.14)

message(STATUS "(EC) Generating mach Makefile")

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_rpn/modules)

include(ec_init)           # Initialise compilers and ec specific functions
ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(mach_VERSION ${PROJECT_VERSION} CACHE INTERNAL "mach version" FORCE)
message(STATUS "(EC) ${PROJECT_NAME} version = ${PROJECT_VERSION}")

enable_language(C Fortran)

ec_build_info()            # Generate build information
include(ec_compiler_presets)

string(TIMESTAMP BUILD_DATE "%Y-%m-%d %H:%M")
if(DEFINED ENV{COMP_ARCH})
  set(EC_COMPILER "$ENV{COMP_ARCH}")
endif()
if(DEFINED ENV{ORDENV_ARCH})
  set(EC_ORDENV_ARCH "$ENV{ORDENV_ARCH}")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/mach_version.in
${CMAKE_CURRENT_SOURCE_DIR}/include/mach_version.inc @ONLY)

#add_subdirectory(src mach)

ec_build_config()                 # Create build configuration script
