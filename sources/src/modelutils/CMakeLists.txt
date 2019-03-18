message(STATUS "Starting generation of modelutils Makefile")

project(MODELUTILS Fortran)

file(GLOB_RECURSE CSRC  *.c)
file(GLOB_RECURSE FTNSRC  *.F90 *.f)

#MESSAGE( STATUS "Fortran compiler: " ${CMAKE_Fortran_COMPILER} )
#MESSAGE( "répertoire des sources: "  ${CMAKE_CURRENT_SOURCE_DIR} )

# if ( NOT "${CMAKE_C_COMPILER}" MATCHES ".*/icc" )
#    if ( "${CMAKE_VERSION}" LESS "3.6" )
#      list(REMOVE_ITEM FTNSRC ${CMAKE_CURRENT_SOURCE_DIR}/ov_ifort/ifort_fpe_handler.F90)
#      list(REMOVE_ITEM FTNSRC ${CMAKE_CURRENT_SOURCE_DIR}/ov_ifort/msg.F90)
#      message( STATUS "Excluding files ifort_fpe_handler.F90 and msg.F90 from sources when compiler = " ${CMAKE_Fortran_COMPILER} )
#    else()
#      list(FILTER FTNSRC EXCLUDE REGEX .*/ov_ifort/.*)
#      list(FILTER FTNSRC EXCLUDE REGEX .*/base/tests/.*)
#      message( STATUS "Excluding directory ov_ifort from sources when compiler = " ${CMAKE_Fortran_COMPILER} )
#    endif()
#  endif()

#MESSAGE( STATUS "FTNSRC: " ${FTNSRC} )

include_directories(AFTER base include)

add_library (modelutils STATIC ${FTNSRC} ${CSRC})
