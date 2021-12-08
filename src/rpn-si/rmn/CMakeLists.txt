message(STATUS "Starting generation of librmn Makefile")

# # gmm requires MPI
# include(FindMPI OPTIONAL)
# set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
# set(CMAKE_C_COMPILER ${MPI_C_COMPILER})

file(GLOB CSRC  */*.c template_utils/whiteboard/*.c)
file(GLOB FTNSRC  */*.F90 */*.F */*.f */*.f90 template_utils/gmm/gmm.f90 template_utils/whiteboard/whiteboard_omp.f90 template_utils/whiteboard/whiteboard_st.f90)

add_library (rmn STATIC ${FTNSRC} ${CSRC})
target_include_directories(rmn PRIVATE
   base
   template_utils
   INTRALIB_INCLUDES
)

if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options( rmn PRIVATE -w )
endif()
