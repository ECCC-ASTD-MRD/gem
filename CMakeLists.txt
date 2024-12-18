cmake_minimum_required(VERSION 3.16)

#----- Append EC specific module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn/modules $ENV{EC_CMAKE_MODULE_PATH})

include(ec_init)           # Include EC specific cmake utils
ec_git_version()           # Get the version from the git repository
ec_parse_manifest()        # Parse MANIFEST file

project(${NAME} DESCRIPTION "${DESCRIPTION}" LANGUAGES C Fortran)
set(PROJECT_VERSION ${VERSION}${STATE})
set(GEM_VERSION ${PROJECT_VERSION} CACHE INTERNAL "GEM version" FORCE)
message(STATUS "(EC) GEM VERSION = ${GEM_VERSION}")

ec_build_info()            # Generate build include file
include(ec_compiler_presets)

# Option for compiling or not with Intel static libraries
option(STATIC_INTEL "Whether to compile with static Intel libraries" OFF)
if ((("$ENV{ECCI_PROCESS}" STREQUAL "staging") OR (STATIC_INTEL)) AND ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel"))
  message(STATUS "(EC) Building with static Intel libraries")
  # Use static MPI libraries: must be set before calling find_package(MPI)
  set(MPI_COMPILER_FLAGS "${MPI_COMPILER_FLAGS} -static_mpi -nostrip")
  set(STATIC_LINK_INTEL_FLAGS "-static-libgcc -static-libstdc++ -qopenmp-link=static")
endif()

find_package(MPI)

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
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qmkl ${STATIC_LINK_INTEL_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qmkl -static-intel -diag-disable 5268 ${STATIC_LINK_INTEL_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -qmkl ${STATIC_LINK_INTEL_FLAGS}")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl ${STATIC_LINK_INTEL_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl -static-intel -diag-disable 5268 ${STATIC_LINK_INTEL_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_INIT "${CMAKE_EXE_LINKER_FLAGS_INIT} -mkl ${STATIC_LINK_INTEL_FLAGS}")
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Warray-temporaries -fdump-core")
  if (EXTRA_CHECKS)
    # in addition to extra checks defined in cmake_rpn compiler presets
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -ffpe-trap=underflow")
  endif()
endif()

if(("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC"))
  find_package(LAPACK)
  find_package(BLAS)
endif()

if((DEFINED ENV{EC_ARCH}) AND ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel"))
  set(FFTW_USE_STATIC_LIBS TRUE)
endif()

find_package(FFTW REQUIRED)
if (WITH_OPENMP)
   set(FFTW FFTW::DoubleOpenMP)
else()
   set(FFTW FFTW::Double)
endif()

# Printing flags
message(STATUS "(EC) CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "(EC) CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")

if (${CMAKE_BUILD_TYPE} MATCHES Debug )
  message(STATUS "(EC) CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "(EC) CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

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
  set(BUILD_ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
  set(EC_ARCH "${BUILD_ARCH}/${EC_COMPILER}")
endif()

if(EXISTS ${CMAKE_SOURCE_DIR}/share/gem-maestro)
  configure_file(${CMAKE_SOURCE_DIR}/share/gem-maestro/version.in ${CMAKE_SOURCE_DIR}/share/gem-maestro/GEM_module/VERSION @ONLY)
  configure_file(${CMAKE_SOURCE_DIR}/share/gem-maestro/version.in ${CMAKE_SOURCE_DIR}/share/gem-maestro/nostage_module/VERSION @ONLY)
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED ENV{GEM_WORK})
    set(WORK_PREFIX "${CMAKE_SOURCE_DIR}/$ENV{GEM_WORK}" CACHE FILEPATH "Working directory prefix" FORCE)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/$ENV{GEM_WORK}" CACHE FILEPATH "Install directory prefix" FORCE)
  else()
    set(WORK_PREFIX "${CMAKE_SOURCE_DIR}/work-${OS}-${COMPILER_SUITE}-${COMPILER_VERSION}" CACHE FILEPATH "Working directory prefix" FORCE)
  endif()
else()
  set(WORK_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE FILEPATH "Working directory prefix" FORCE)
endif()
message(STATUS "(EC) Setting working directory prefix to ${WORK_PREFIX}")

add_subdirectory(${CMAKE_SOURCE_DIR}/src src)

# list of binaries depends on built sources
list(APPEND BINARIES maingemdm)
if((NOT DEFINED ENV{EC_ARCH}) OR (NOT WITH_SYSTEM_RPN))
  list(APPEND BINARIES cclargs_lite r.fstinfo r.filetype editfst fstcomp fststat pgsm voir )
  list(APPEND LINKS ${CMAKE_BINARY_DIR}/src/fst-tools/src/editfst/editfst)
  list(APPEND LINKS ${CMAKE_BINARY_DIR}/src/fst-tools/src/fstcomp/fstcomp)
  list(APPEND LINKS ${CMAKE_BINARY_DIR}/src/fst-tools/src/fststat/fststat)
  list(APPEND LINKS ${CMAKE_BINARY_DIR}/src/fst-tools/src/pgsm/pgsm)
  list(APPEND LINKS ${CMAKE_BINARY_DIR}/src/fst-tools/src/voir/voir)
  install(FILES ${LINKS} DESTINATION ${WORK_PREFIX}/bin)
endif()
list(APPEND BINARIES yy2global flipit yydecode yyencode test_integrals)
if ("${phy_LIBRARIES}" STREQUAL "rpnphy" )
  list(APPEND BINARIES feseri prphynml )
endif()
list(APPEND BINARIES gemgrid checkdmpart prgemnml toc2nml gem_monitor_output gem_monitor_end)

install(TARGETS ${BINARIES} RUNTIME DESTINATION ${WORK_PREFIX}/bin)

# install only a subset of targets
add_custom_target(work
    COMMAND /bin/mkdir -p ${WORK_PREFIX}
    COMMAND /bin/mkdir -p ${WORK_PREFIX}/bin
    COMMAND /bin/ln -sf ${CMAKE_SOURCE_DIR}/share/configurations ${WORK_PREFIX}
    COMMAND cd ${CMAKE_BINARY_DIR} && make install/local
  )

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
