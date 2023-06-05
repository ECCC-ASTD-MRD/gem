cmake_minimum_required(VERSION 3.14)

message(STATUS "(EC) Generating rpnphy Makefile")

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_rpn/modules)

include(ec_init)           # Initialise compilers and ec specific functions
ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(rpnphy_VERSION ${PROJECT_VERSION} CACHE INTERNAL "rpnphy version" FORCE) # Needed for cascaded version identification
message(STATUS "(EC) ${PROJECT_NAME} version = ${PROJECT_VERSION}")

enable_language(C Fortran)

ec_build_info()            # Generate build information
include(ec_compiler_presets)

# To build without OpenMP, you have to add -DWITH_OPENMP=FALSE to the cmake command line
set(WITH_OPENMP TRUE CACHE BOOL "Control whether to use OpenMP")
include(ec_openmp)

if (NOT rmn_FOUND)
  find_package(rmn ${rmn_REQ_VERSION} COMPONENTS static ompi REQUIRED)
endif()

if (NOT rpncomm_FOUND)
  find_package(rpncomm ${rpncomm_REQ_VERSION} REQUIRED)
endif()

if (NOT tdpack_FOUND)
  find_package(tdpack ${tdpack_REQ_VERSION} COMPONENTS static REQUIRED)
endif()

if (NOT vgrid_FOUND)
  find_package(vgrid ${vgrid_REQ_VERSION} COMPONENTS static REQUIRED)
endif()

# Adding specific flags for rpnphy
if (("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel") AND NOT ("${CMAKE_SYSTEM_NAME}" STREQUAL "CrayLinuxEnvironment"))
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 2021)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qmkl")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -static-intel -diag-disable 5268")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -qmkl")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -static-intel -diag-disable 5268")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -mkl")
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  find_package(LAPACK)
  find_package(BLAS)
endif()

string(TIMESTAMP BUILD_DATE "%Y-%m-%d %H:%M")
if(DEFINED ENV{COMP_ARCH})
  set(EC_COMPILER "$ENV{COMP_ARCH}")
endif()
if(DEFINED ENV{ORDENV_ARCH})
  set(EC_ORDENV_ARCH "$ENV{ORDENV_ARCH}")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/rpnphy_version.in ${CMAKE_CURRENT_SOURCE_DIR}/include/rpnphy_version.inc @ONLY)

add_subdirectory(src/main rpnphy-main)
add_subdirectory(src rpnphy)

ec_build_config()                 # Create build configuration script
