message(STATUS "(EC) Generating rpnphy libraries Makefile")

find_package(MPI)
if (MPI_FOUND)
   set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif()

file(GLOB F_FILES main/*.F90 api/*.F90 base/*.F90 series/*.F90 surface/*.F90 surface/from_surfex/*.F90 utils/*.F90 ${MACH_FILES})
list(FILTER F_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/src/rpnphy/main/main*")
if(WITH_NEMO)
  file(GLOB CPLOCN_SOURCES ${CPLOCN_GLOB})
  list(APPEND F_FILES ${CPLOCN_SOURCES})
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
set(PROJECT_INCLUDE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(rpnphy STATIC ${F_FILES})

target_compile_definitions(rpnphy PUBLIC ECCCGEM CROCUS_EXT)
if (WITH_MACH)
   target_compile_definitions(rpnphy PUBLIC HAVE_MACH)
   message(STATUS "(EC) rpnphy compile definition: HAVE_MACH")
   target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${MACH_INCLUDE_DIRECTORY}>)
endif()
if(WITH_NEMO)
  target_compile_definitions(rpnphy PUBLIC HAVE_NEMO)
  message(STATUS "(EC) rpnphy compile definition: HAVE_NEMO")
endif()

target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${PROJECT_INCLUDE_DIRECTORY}>)
target_include_directories(rpnphy PUBLIC $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

target_include_directories(rpnphy PRIVATE
  utils
  base
  )

target_link_libraries(rpnphy PUBLIC
  OpenMP::OpenMP_Fortran
  OpenMP::OpenMP_C
  vgrid::vgrid
  tdpack::tdpack
  modelutils
  )

if(WITH_NEMO)
  string(REPLACE ":" ";" LIBRARY_PATH_AS_LIST $ENV{LIBRARY_PATH})
  list(APPEND CMAKE_LIBRARY_PATH ${LIBRARY_PATH_AS_LIST})
  find_library(LIBCPLGEMNEMO cplgemnemo REQUIRED)
  target_link_libraries(rpnphy PUBLIC ${LIBCPLGEMNEMO})
endif()

set_target_properties(rpnphy PROPERTIES VERSION ${PROJECT_VERSION})
add_dependencies(rpnphy ${PROJECT_NAME}_build_info)


install(TARGETS rpnphy
    LIBRARY       DESTINATION lib
    ARCHIVE       DESTINATION lib)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
