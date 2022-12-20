cmake_minimum_required(VERSION 3.10)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH $ENV{EC_CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_rpn/modules)

include(ec_init)           # Include EC specific cmake utils
ec_git_version()           # Get the version from the git repository
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}")
set(PROJECT_VERSION ${VERSION}${STATE})
set(GEM_VERSION ${PROJECT_VERSION} CACHE INTERNAL "GEM version" FORCE)
message(STATUS "(EC) GEM VERSION = ${GEM_VERSION}")

enable_language(Fortran)
enable_language(C)

ec_build_info()            # Generate build include file
include(ec_compiler_presets)

# To build without OpenMP, you have to add -DWITH_OPENMP=FALSE to the cmake command line
set(WITH_OPENMP TRUE CACHE BOOL "Control whether to use OpenMP")
include(ec_openmp)

option(WITH_OPENACC "Build with OpenACC support" OFF)
if (WITH_OPENACC)
   find_package(OpenACC)
   if (OPENACC_NOTFOUND)
      message(FATAL_ERROR "(EC) OpenACC was requested, but was not found!")
   endif()
   if(OpenACC_C_FLAGS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS} ${OpenACC_extra_FLAGS}")
   endif()
   if(OpenACC_Fortran_FLAGS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_Fortran_FLAGS}")
   elseif(OpenACC_C_FLAGS)
      # Fallback on C flags and hope it also works for Fortran
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_C_FLAGS} ${OpenACC_extra_FLAGS}")
   endif()
else()
   message(STATUS "(EC) Building WITHOUT OpenACC")
endif()

option(EXTRA_CHECKS "Compile with extra checks to detect things such as overunned bounds and stack trashing" OFF)
# Adding specific flags for GEM
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
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Warray-temporaries -fdump-core")
  if (EXTRA_CHECKS)
    # in addition to extra checks defined in cmake_rpn compiler presets
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=underflow -fimplicit-none")
  endif()
endif()

if(("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC"))
  find_package(LAPACK)
  find_package(BLAS)
endif()

# Printing flags
message(STATUS "(EC) CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "(EC) CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")

if (${CMAKE_BUILD_TYPE} MATCHES Debug )
  message(STATUS "(EC) CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "(EC) CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

# Figure out the correct arch-dependent include path
set(BUILD_ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "(EC) Build architecture: ${BUILD_ARCH}")
message(STATUS "(EC) Sources directory is: ${CMAKE_SOURCE_DIR}")
message(STATUS "(EC) Build directory is: ${CMAKE_BINARY_DIR}")

# Instead of adding -I. to the FLAGS, we can let cmake add the current folder automatically
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# To see the make commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

# get name and version of operating system
execute_process(COMMAND sh "-c" "${CMAKE_SOURCE_DIR}/utils/os.sh" OUTPUT_VARIABLE OS)
message(STATUS "(EC) Operating system is: ${OS}")

# get name and version of compiler
if(DEFINED ENV{COMPILER_VERSION})
   set(COMPILER_VERSION $ENV{COMPILER_VERSION})
else()
  execute_process(COMMAND sh "-c" "${CMAKE_SOURCE_DIR}/utils/compiler.sh ${COMPILER_SUITE}" OUTPUT_VARIABLE COMPILER_VERSION)
endif()
message(STATUS "(EC) Compiler version: ${COMPILER_VERSION}")

# Some variables for cmake-configured files
string(TIMESTAMP BUILD_DATE "%Y-%m-%d %H:%M")
if(DEFINED ENV{COMP_ARCH})
  set(EC_COMPILER "$ENV{COMP_ARCH}")
else()
  set(EC_COMPILER "${COMPILER_SUITE}-${COMPILER_VERSION}")
endif()
if(DEFINED ENV{EC_PLAT})
  set(EC_PLAT "$ENV{EC_PLAT}")
else()
  set(EC_PLAT "${OS}")
endif()
if(DEFINED ENV{ORDENV_ARCH})
  set(EC_ORDENV_ARCH "$ENV{ORDENV_ARCH}")
endif()
if(DEFINED ENV{EC_ARCH})
  set(EC_ARCH "$ENV{EC_ARCH}")
else()
  set(EC_ARCH "${BUILD_ARCH}/${EC_COMPILER}")
endif()

set(WORK_PREFIX ${ROOT} CACHE FILEPATH "Working environment directory")
if(DEFINED ENV{GEM_WORK})
    set(WORK_PREFIX "${CMAKE_SOURCE_DIR}/$ENV{GEM_WORK}" CACHE FILEPATH "Working directory prefix" FORCE)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/$ENV{GEM_WORK}" CACHE FILEPATH "Install directory prefix" FORCE)
else()
    set(WORK_PREFIX "${CMAKE_SOURCE_DIR}/work-${OS}-${COMPILER_SUITE}-${COMPILER_VERSION}" CACHE FILEPATH "Working directory prefix" FORCE)
endif()
message(STATUS "(EC) Setting working directory prefix to ${WORK_PREFIX}")

add_subdirectory(${CMAKE_SOURCE_DIR}/src src)

# list of binaries depends on built sources
list(APPEND BINARIES maingemdm cclargs_lite fstinfo r.filetype)
if (NOT DEFINED ENV{EC_ARCH})
  list(APPEND BINARIES editfst fstcomp fststat pgsm voir )
endif()
list(APPEND BINARIES yy2global flipit yydecode yyencode)
if ("${phy_SRC}" STREQUAL "rpnphy" )
  list(APPEND BINARIES feseri prphynml )
endif()
list(APPEND BINARIES gemgrid checkdmpart prgemnml toc2nml gem_monitor_output gem_monitor_end )

install(TARGETS ${BINARIES} RUNTIME DESTINATION ${WORK_PREFIX}/bin)

add_custom_target(work
    COMMAND /bin/mkdir -p ${WORK_PREFIX}
    COMMAND /bin/mkdir -p ${WORK_PREFIX}/bin
    COMMAND /bin/ln -sf ${CMAKE_SOURCE_DIR}/configurations ${WORK_PREFIX}
    COMMAND /bin/ln -sf ${WORK_PREFIX}/bin/fstinfo ${WORK_PREFIX}/bin/r.fstinfo
    COMMAND cd ${CMAKE_BINARY_DIR} && make install/local
  )

ec_package_name()      # Define package prefix
ec_build_config()      # Create build configuration script
ec_prepare_ssm()       # Prepare ssm packaging files

set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "ECCC")
set(CPACK_PACKAGE_CONTACT "${MAINTAINER}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_OUTPUT_FILE_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/package")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${NAME}_${GEM_VERSION}")
include(CPack)
