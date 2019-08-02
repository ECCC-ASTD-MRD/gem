# - Find FFTW
# Find the fftw includes and library
#
#  FFTW_INCLUDE_DIR - where to find fftw3.h, fftw3.f03, etc.
#  FFTW_LIBRARIES   - List of libraries when using fftw.
#  FFTW_FOUND       - True if fftw found.

IF (FFTW_INCLUDE_DIR)
  # Already in cache, be silent
  SET(FFTW_FIND_QUIETLY TRUE)
ENDIF()

FIND_PATH(FFTW_INCLUDE_DIR
  NAMES fftw3.h
  PATHS /usr/include
        ENV CPATH
)
INCLUDE_DIRECTORIES(${FFTW_INCLUDE_DIR})

# TODO: This will need to be modified to work with the system
# provided fftw on p9gpu-*
# On RHEL, librairies have a version suffix

find_library(FFTW_LIB_PATH
   NAMES fftw3
   PATHS /usr/lib64
         /usr/lib
         ENV LD_LIBRARY_PATH
)
find_library(FFTW_OMP_LIB_PATH
   NAMES fftw3_omp
   PATHS /usr/lib64
         /usr/lib
         ENV LD_LIBRARY_PATH
)
list(APPEND FFTW_LIBRARIES ${FFTW_LIB_PATH} ${FFTW_OMP_LIB_PATH})

if (FFTW_INCLUDE_DIR AND FFTW_LIB_PATH AND FFTW_OMP_LIB_PATH)
   set(FFTW_FOUND TRUE)
endif()

message(STATUS "FFTW_INCLUDE_DIR : ${FFTW_INCLUDE_DIR}")
message(STATUS "FFTW_LIBRARIES : ${FFTW_LIBRARIES}")
