cmake_minimum_required(VERSION 3.14)

message(STATUS "(EC) Generating modelutils Makefile")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn)
include(ec_init)           # Initialise compilers and ec specific functions
ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(modelutils_VERSION ${PROJECT_VERSION} CACHE INTERNAL "modelutils version" FORCE) # Needed for cascaded version identification
message(STATUS "(EC) ${PROJECT_NAME} version = ${PROJECT_VERSION}")

enable_language(C Fortran)

ec_build_info()            # Generate build information
include(ec_compiler_presets)
include(ec_openmp)

# necessary for rmn 19.7.0
if (OPENMP_FOUND)
  add_compile_definitions(WITH_OpenMP)
else()
  add_compile_definitions(WITHOUT_OpenMP)
endif()

if (NOT RMN_FOUND)
  find_package(RMN ${RMN_REQ_VERSION} REQUIRED)
endif()

if (NOT rpncomm_FOUND)
   find_package(rpncomm ${rpncomm_REQ_VERSION} CONFIG)
endif()

if (NOT tdpack_FOUND)
  find_package(tdpack ${tdpack_REQ_VERSION} CONFIG)
endif()

if (NOT VGRID_FOUND)
  find_package(VGRID ${VGRID_REQ_VERSION})
endif()

  # Adding specific flags for modelutils
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
  set(LAPACK_LIBRARIES "lapack")
  set(BLAS_LIBRARIES "blas")
endif()

string(TIMESTAMP BUILD_DATE "%Y-%m-%d %H:%M")
if(DEFINED ENV{COMP_ARCH})
  set(EC_COMPILER "$ENV{COMP_ARCH}")
endif()
if(DEFINED ENV{ORDENV_ARCH})
  set(EC_ORDENV_ARCH "$ENV{ORDENV_ARCH}")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/modelutils_version.in ${CMAKE_CURRENT_SOURCE_DIR}/include/modelutils_version.inc @ONLY)

add_subdirectory(src/main modelutils-main)
add_subdirectory(src modelutils)


#----- Generate the config file for the project to be usable via cmake's find_package command
set(INCLUDE_INSTALL_DIR include)
set(LIB_INSTALL_DIR     lib)
set(CONFIG_INSTALL_DIR  "${LIB_INSTALL_DIR}/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${CONFIG_INSTALL_DIR}"
    PATH_VARS           INCLUDE_INSTALL_DIR LIB_INSTALL_DIR
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion
)
install(FILES   "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${CONFIG_INSTALL_DIR}"
)

#----- Packaging
ec_package_name()                 # Define package prefix  
ec_build_config()                 # Create build configuration script
ec_prepare_ssm()                  # Prepare ssm packaging files

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "ECCC")
set(CPACK_PACKAGE_CONTACT "${MAINTAINER}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/package")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${NAME}_${PROJECT_VERSION}")
include(CPack)

if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${PACKAGE_NAME}" CACHE FILEPATH "CMake Installation prefix" FORCE)
endif()
