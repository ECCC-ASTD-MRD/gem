cmake_minimum_required(VERSION 3.14)

message(STATUS "(EC) Generating modelutils Makefile")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn)
include(ec_init)           # Initialise compilers and ec specific functions
#ec_git_version()           # Get version from git state
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(MODELUTILS_VERSION ${PROJECT_VERSION} CACHE INTERNAL "Modelutils version" FORCE) # Needed for cascaded version identification

enable_language(C Fortran)

#ec_build_info()            # Generate build information
ec_install_prefix()        # Define install prefix  

include(ec_compiler_presets)
include(ec_openmp)           # Enable OpenMP (-DWITH_OPENMP=TRUE)

if (NOT RMN_FOUND)
  find_package(RMN ${RMN_REQ_VERSION})
endif()

if (NOT VGRID_FOUND)
  find_package(VGRID ${VGRID_REQ_VERSION})
endif()

# Adding specific flags for modelutils
if (("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel") AND NOT ("${CMAKE_SYSTEM_NAME}" STREQUAL "CrayLinuxEnvironment"))
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -static-intel -diag-disable 5268")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -mkl")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(LAPACK_LIBRARIES "lapack")
  set(BLAS_LIBRARIES "blas")
endif()

add_subdirectory(main modelutils-main)

include(FindMPI)
set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB C_FILES  */*.c)
file(GLOB F_FILES  main/flipit.F90 main/yy2global.F90 main/yydecode.F90 main/yyencode.F90 base/*F90 driver/*.F90 utils/*.F90 tdpack/*F90 stubs/*F90 tmg_stubs/*F90 ov_ifort/*F90)

add_library (modelutils STATIC ${F_FILES} ${C_FILES})
#add_dependencies(modelutils ${PROJECT_NAME}_build_info)
target_include_directories(modelutils PRIVATE base include)

add_library (modelutils_ov_ifort STATIC ov_ifort/ifort_fpe_handler.F90)
add_dependencies(modelutils_ov_ifort modelutils)
#add_dependencies(modelutils_ov_ifort ${PROJECT_NAME}_build_info)
target_include_directories(modelutils_ov_ifort PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

add_library (modelutils_stubs STATIC stubs/drv_dyn_itf_stubs.F90)
add_dependencies(modelutils_stubs modelutils)
#add_dependencies(modelutils_stubs ${PROJECT_NAME}_build_info)

add_library (modelutils_tmg_stubs STATIC tmg_stubs/tmg_stubs.F90)
#add_dependencies(modelutils_tmg_stubs ${PROJECT_NAME}_build_info)

ec_build_config()      # Create build configuration script
ec_prepare_ssm()       # Prepare ssm packaging files

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "ECCC")
set(CPACK_PACKAGE_CONTACT "${MAINTAINER}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
#set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.org")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/package")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_FILE_NAME "${CMAKE_INSTALL_PREFIX}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_INSTALL_PREFIX}")
include(CPack)

install(TARGETS modelutils modelutils)
install(TARGETS modelutils_ov_ifort modelutils_ov_ifort)
install(TARGETS modelutils_stubs modelutils_stubs)
install(TARGETS modelutils_tmg_stubs modelutils_tmg_stubs)
install(PROGRAMS ${CMAKE_BINARY_DIR}/modelutils-config DESTINATION ${CMAKE_INSTALL_PREFIX})
