cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_rpn)
message(STATUS "CMAKE_MODULE_PATH=" ${CMAKE_SOURCE_DIR}/cmake_rpn)

include(ec_init)           # Include EC specific cmake utils
ec_git_version()           # Get the version from the git repository

project(GEM DESCRIPTION "GEM - Forecasting and data assimilation system")
#ec_build_info()            # Generate build include file (optional)

enable_language(Fortran)
enable_language(C)

include(ec_compiler_presets)

# CMAKE_BUILD_TYPE can be one of Debug, Release, RelWithDebInfo, MinSizeRel
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "CMAKE_BUILD_TYPE can be one of Debug, Release, RelWithDebInfo, MinSizeRel")
  message(STATUS "No build type selected, default to RelWithDebInfo")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()
message(STATUS "Configuring for ${CMAKE_BUILD_TYPE} build type")

# To build without OpenMP, you have to add -DWITH_OPENMP=FALSE to the cmake command line
set(WITH_OPENMP TRUE CACHE BOOL "Control whether to use OpenMP")
include(OpenMP)

option(WITH_OPENACC "Build with OpenACC support" OFF)
if (WITH_OPENACC)
   find_package(OpenACC)
   if (OPENACC_NOTFOUND)
      message(FATAL_ERROR "OpenACC was requested, but was not found!")
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
   message(STATUS "Building WITHOUT OpenACC")
endif()

option(EXTRA_CHECKS "Compile with extra checks to detect things such as overunned bounds and stack trashing" OFF)

if (${CMAKE_BUILD_TYPE} MATCHES Debug )
  message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
endif()

# Figure out the correct arch-dependent include path
set(BUILD_ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Build architecture: ${BUILD_ARCH}")
#message(STATUS "Basic CMAKE_C_FLAGS (cmake_rpn): ${CMAKE_C_FLAGS}")
#message(STATUS "Basic CMAKE_Fortran_FLAGS (cmake_rpn): ${CMAKE_Fortran_FLAGS}")

message(STATUS "Sources directory is: ${CMAKE_SOURCE_DIR}")
message(STATUS "Build directory is: ${CMAKE_BINARY_DIR}")

file(READ "${CMAKE_SOURCE_DIR}/VERSION" GEM_VERSION)
string(REGEX REPLACE "\n$" "" GEM_VERSION "${GEM_VERSION}")
message(STATUS "GEM_VERSION: ${GEM_VERSION}")

# Include an architecture dependent file that contains include files for a
# particular architecture and compiler, using the default value of the
# compiler or the argument given to cmake command.
message(STATUS "arch-dependent include directory: ${CMAKE_SOURCE_DIR}/src/rpn-si/archdep/include/${BUILD_ARCH}-${COMPILER_SUITE}")
include_directories(
   ${CMAKE_SOURCE_DIR}/src/rpn-si/archdep/include/${BUILD_ARCH}-${COMPILER_SUITE}
   ${CMAKE_SOURCE_DIR}/src/rpn-si/mpi-extras
   ${CMAKE_SOURCE_DIR}/src/rpn-si/rpncomm
   ${CMAKE_SOURCE_DIR}/src/rpn-si/rmn/PUBLIC_INCLUDES
   ${CMAKE_SOURCE_DIR}/src/rpn-si/RPN/compiler_includes
)

# Instead of adding -I. to the FLAGS, we can let cmake add the current folder automatically
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# To see the make commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

# get name and version of operating system
execute_process(COMMAND sh "-c" "${CMAKE_SOURCE_DIR}/project/os.sh" OUTPUT_VARIABLE OS)
message(STATUS "Operating system is: ${OS}")

# get name and version of compiler
if(DEFINED ENV{COMPILER_VERSION})
   set(COMPILER_VERSION $ENV{COMPILER_VERSION})
else()
  execute_process(COMMAND sh "-c" "${CMAKE_SOURCE_DIR}/project/compiler.sh ${COMPILER_SUITE}" OUTPUT_VARIABLE COMPILER_VERSION)
endif()
message(STATUS "COMPILER_VERSION: ${COMPILER_VERSION}")

if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local")
  if(DEFINED ENV{GEM_WORK})
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/$ENV{GEM_WORK}" CACHE FILEPATH "CMake Installation prefix" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/work-${OS}-${COMPILER_SUITE}-${COMPILER_VERSION}" CACHE FILEPATH "CMake Installation prefix" FORCE)
  endif()
message(STATUS "Setting CMAKE_INSTALL_PREFIX to ${CMAKE_INSTALL_PREFIX}")
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/src src)

add_custom_target(work
    /bin/mkdir -p ${CMAKE_INSTALL_PREFIX}/bin
    COMMAND /bin/cp -avf
        $<TARGET_FILE:cclargs_lite>
        $<TARGET_FILE:editfst>
        $<TARGET_FILE:flipit>
        $<TARGET_FILE:fstcomp>
        $<TARGET_FILE:fststat>
        $<TARGET_FILE:maingemdm>
        $<TARGET_FILE:gemgrid>
        $<TARGET_FILE:yy2global>
        $<TARGET_FILE:yydecode>
        $<TARGET_FILE:yyencode>
        $<TARGET_FILE:feseri>
        $<TARGET_FILE:prphynml>
        $<TARGET_FILE:checkdmpart>
        $<TARGET_FILE:prgemnml>
        $<TARGET_FILE:toc2nml>
        $<TARGET_FILE:gem_monitor_output>
        $<TARGET_FILE:gem_monitor_end>
        $<TARGET_FILE:pgsm>
        $<TARGET_FILE:r.fstinfo>
        $<TARGET_FILE:r.filetype>
        $<TARGET_FILE:voir>
        ${CMAKE_INSTALL_PREFIX}/bin
    COMMAND /bin/ln -sf ${CMAKE_SOURCE_DIR}/configurations ${CMAKE_INSTALL_PREFIX}
    COMMAND /bin/cp -avf ${CMAKE_SOURCE_DIR}/doc/README.run ${CMAKE_INSTALL_PREFIX}

    DEPENDS
        cclargs_lite
        flipit
        maingemdm
        gemgrid
        yy2global
	yyencode
	yydecode
        feseri
        prphynml
        checkdmpart
        prgemnml
        toc2nml
        gem_monitor_output
        gem_monitor_end
        r.fstinfo
        r.filetype
        fstcomp
        fststat
        editfst
        voir
)
