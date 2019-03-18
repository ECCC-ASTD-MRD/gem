message(STATUS "Starting generation of librmn Makefile")

project(RMN Fortran)

file(GLOB CSRC  */*.c template_utils/whiteboard/*.c)
file(GLOB FTNSRC  */*.F90 */*.F */*.f */*.f90 template_utils/gmm/gmm.f90 template_utils/whiteboard/whiteboard_omp.f90 template_utils/whiteboard/whiteboard_st.f90)

include_directories(AFTER base template_utils INTRALIB_INCLUDES)

add_library (rmn STATIC ${FTNSRC} ${CSRC})
#add_library (rmn STATIC ${OUTPUT})
