cmake_minimum_required(VERSION 3.10)

project(GEM DESCRIPTION "GEM - Forecasting and data assimilation system" LANGUAGES NONE)

cmake_policy(SET CMP0066 NEW)
cmake_policy(SET CMP0056 NEW)

# set the compiler suite by uncommenting one of the following 3 lines:
set(COMPILER_SUITE intel CACHE STRING compiler)
#set(COMPILER_SUITE gnu CACHE STRING compiler)
#set(COMPILER_SUITE pgi CACHE STRING compiler)
message(STATUS "Compiler suite to be used is: ${COMPILER_SUITE}")

# CMAKE_BUILD_TYPE can be one of Debug, Release, RelWithDebInfo, MinSizeRel
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "CMAKE_BUILD_TYPE can be one of Debug, Release, RelWithDebInfo, MinSizeRel")
  message(STATUS "No build type selected, default to RelWithDebInfo")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()
message(STATUS "Configuring for ${CMAKE_BUILD_TYPE} build type")

# Include a file that contains settings for a particular architecture and
# compiler, using the default value of the compiler defined above or the
# argument given to cmake command.
set(BUILD_ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

if(DEFINED ENV{ORDENV_PLAT})
   message(STATUS "include ${CMAKE_SOURCE_DIR}/project/$ENV{ORDENV_PLAT}-${COMPILER_SUITE}.cmake")
   include("${CMAKE_SOURCE_DIR}/project/$ENV{ORDENV_PLAT}-${COMPILER_SUITE}.cmake")
else()
  message(STATUS "include ${CMAKE_SOURCE_DIR}/project/${BUILD_ARCH}-${COMPILER_SUITE}.cmake")
  include("${CMAKE_SOURCE_DIR}/project/${BUILD_ARCH}-${COMPILER_SUITE}.cmake")
endif()

option(EXTRA_CHECKS "Compile with extra checks to detect things such as overunned bounds and stack trashing" OFF)

# To build without OpenMP, you have to add -DWITH_OPENMP=0 to the cmake command line
option(WITH_OPENMP "Build with OpenMP support" ON)
if (WITH_OPENMP)
   find_package(OpenMP)
   if (NOT OPENMP_FOUND)
      message(FATAL_ERROR "OpenMP was requested, but was not found!")
   endif()
   if(OpenMP_C_FLAGS)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   endif()
   if(OpenMP_Fortran_FLAGS)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
   else()
      # If we can't find Fortran flags, fallback on C flags
      message(WARNING "OpenMP couldn't find Fortran OpenMP Flags.  Falling back on C Flags.")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
   endif()
else()
   message(STATUS "Building WITHOUT OpenMP")
endif()

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

message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")
if (${CMAKE_BUILD_TYPE} MATCHES Debug )
  message(STATUS "CMAKE_C_FLAGS_DEBUG=${CMAKE_C_FLAGS_DEBUG}")
  message(STATUS "CMAKE_Fortran_FLAGS_DEBUG=${CMAKE_Fortran_FLAGS_DEBUG}")
endif()


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

if(DEFINED ENV{GEM_WORK})
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/$ENV{GEM_WORK}" CACHE FILEPATH "CMake Installation prefix" FORCE)
else()
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/work-${OS}-${COMPILER_SUITE}-${COMPILER_VERSION}" CACHE FILEPATH "CMake Installation prefix" FORCE)
endif()
message(STATUS "Setting CMAKE_INSTALL_PREFIX to ${CMAKE_INSTALL_PREFIX}")

add_subdirectory(${CMAKE_SOURCE_DIR}/src src)

add_custom_target(work
    /bin/mkdir -p ${CMAKE_INSTALL_PREFIX}/bin
    COMMAND /bin/cp -avf
        $<TARGET_FILE:cclargs_lite>
        $<TARGET_FILE:flipit>
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
        $<TARGET_FILE:r.fstinfo>
        $<TARGET_FILE:r.filetype>
        ${CMAKE_INSTALL_PREFIX}/bin
    COMMAND /bin/ln -sf ${CMAKE_SOURCE_DIR}/configurations ${CMAKE_INSTALL_PREFIX}

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
)
