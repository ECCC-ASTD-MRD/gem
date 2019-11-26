cmake_minimum_required(VERSION 2.8.12)

project(GEM NONE)

# set the compiler suite by uncommenting one of the following 3 lines:
set(COMPILER gnu CACHE STRING compiler)
#set(COMPILER intel CACHE compiler)
#set(COMPILER pgi CACHE compiler)
# which will set the compiler and you can just type:
message(STATUS "Compiler to be used is: ${COMPILER}")

# CMAKE_BUILD_TYPE can be one of Debug, Release, RelWithDebInfo, MinSizeRel
set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)")
message(STATUS "Configuring for ${CMAKE_BUILD_TYPE}")

# Include an architecture dependent file that contains settings for a
# particular architecture and compiler, using the default value of the
# compiler defined above or the argument given to cmake command.
set(ARCH "${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")

message(STATUS "include ${ARCH}-${COMPILER}.cmake")
include("${ARCH}-${COMPILER}.cmake")


# TODO: This isn't used presently.  Check to see if we need it
# set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

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
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
   endif()
else()
   message(STATUS "Building WITHOUT OpenMP")
endif()

message(STATUS "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}")

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

# get name and version of operating system
execute_process(COMMAND sh "-c" "${CMAKE_CURRENT_SOURCE_DIR}/os.sh" OUTPUT_VARIABLE OS)
message(STATUS "Operating system is: ${OS}")

# get name and version of compiler
execute_process(COMMAND sh "-c" "${CMAKE_CURRENT_SOURCE_DIR}/compiler.sh ${COMPILER}" OUTPUT_VARIABLE COMP_VERSION)
message(STATUS "Compiler version is: ${COMPILER} ${COMP_VERSION}")

# Enable distclean target
include("distclean.cmake")

set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;f90;F90 CACHE STRING "Fortran extension files")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/build/modules CACHE STRING "Where Fortran modules go")

message(STATUS "Sources directory is: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Build directory is: ${CMAKE_CURRENT_BINARY_DIR}")

get_filename_component(ROOT ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY CACHE)
message(STATUS "ROOT=${ROOT}")

# Include an architecture dependent file that contains include files for a
# particular architecture and compiler, using the default value of the
# compiler defined above or the argument given to cmake command.
message(STATUS "arch-dependent include directory: ${ROOT}/src/archdep/include/${ARCH}-${COMPILER}")

include_directories(
   ${ROOT}/src/archdep/include/${ARCH}-${COMPILER}
   ${CMAKE_CURRENT_BINARY_DIR}/build/modules
   ${ROOT}/src/rpncomm
   ${ROOT}/src/rmn/PUBLIC_INCLUDES
   ${ROOT}/src/RPN/compiler_includes
)

# Instead of adding -I. to the FLAGS, we can let cmake add the current folder automatically
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Too see the make commands
#set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory(${ROOT}/src build/src)

install(
   DIRECTORY ${ROOT}/scripts/scripts
   DESTINATION ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION}
   USE_SOURCE_PERMISSIONS
)
install(
   FILES ${ROOT}/scripts/runprep
   DESTINATION ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION}
   PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
install(
   FILES ${ROOT}/scripts/runmod
   DESTINATION ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION}
   PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
install(
   FILES ${ROOT}/scripts/mybatchjob.pbs
      DESTINATION ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION}
      PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)
install(
   FILES ${ROOT}/doc/README.run
   DESTINATION ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION}
   PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)
install(
   CODE "execute_process(COMMAND
      /bin/ln -sf ${ROOT}/work/configurations ${ROOT}/work/work-${OS}-${COMPILER}-${COMP_VERSION} ) "
)
